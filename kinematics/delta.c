/*
  delta.c - delta kinematics implementation

  Part of grblHAL

  Copyright (c) 2023 Terje Io
  Transforms derived from mzavatsky at Trossen Robotics
    https://hypertriangle.com/~alex/delta-robot-tutorial/
  get_cuboid_envelope() derived from javascript code in
    https://www.marginallyclever.com/other/samples/fk-ik-test.html

  Grbl is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Grbl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Grbl.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../grbl.h"

#if DELTA_ROBOT

#include <math.h>
#include <string.h>

#include "../hal.h"
#include "../settings.h"
#include "../nvs_buffer.h"
#include "../planner.h"
#include "../motion_control.h"
#include "../protocol.h"
#include "../kinematics.h"

#define A_MOTOR X_AXIS
#define B_MOTOR Y_AXIS
#define C_MOTOR Z_AXIS

// Define step segment generator state flags.
typedef union {
    uint8_t mask;
    struct {
        uint8_t home_to_cuboid_top :1,
                limit_to_cuboid    :1,
                unassigned         :6;
    };
} delta_settings_flags_t;

typedef struct {
    float rf;   // bicep length
    float re;   // forearm length
    float f;    // base radius
    float e;    // end effector radius
    float b;    // base to floor
    float sl;   // max segment length TODO: replace with calculated value?
    float home_angle;
    float max_angle;
    delta_settings_flags_t flags;
} delta_settings_t;

typedef struct {
    float home_z;
    float home_pulloff;
    float resolution;
    float envelope_midz;
    float min_angle[3];
    float max_angle[3];
    coord_data_t last_pos;
    float t;
    float y1;
    float y1_sqr;
    float fe;
    float rf_sqr;   // bicep length ^ 2
    float re_sqr;   // forearm length ^ 2
    delta_settings_t cfg;
} delta_properties_t;

static bool jog_cancel = false;
static delta_settings_t delta_settings;
static delta_properties_t machine = {0};
static homing_mode_t homing_mode;
static on_report_options_ptr on_report_options;
static on_homing_completed_ptr on_homing_completed;
static on_report_command_help_ptr on_report_command_help;
static on_set_axis_setting_unit_ptr on_set_axis_setting_unit;
static on_setting_get_description_ptr on_setting_get_description;
static jog_limits_ptr apply_jog_limits;
static nvs_address_t nvs_address;

// inverse kinematics
// helper functions, calculates angle position[X_AXIS] (for YZ-pane)
int delta_calcAngleYZ (float x0, float y0, float z0, float *theta)
{
    y0 -= machine.fe;    // shift center to edge
    // z = a + b*y
    float a = (x0 * x0 + y0 * y0 + z0 * z0 + machine.rf_sqr - machine.re_sqr - machine.y1_sqr) / (2.0f * z0);
    float b = (machine.y1 - y0) / z0;
    // discriminant
    float d = -(a + b * machine.y1) * (a + b * machine.y1) + machine.cfg.rf * (b * b * machine.cfg.rf + machine.cfg.rf);
    if (d < 0.0f)
        return -1; // non-existing point

    float yj = (machine.y1 - a * b - sqrt(d)) / (b * b + 1); // choosing outer point
    float zj = a + b * yj;
 //   *theta = 180.0f * atanf(-zj / (y1 - yj)) / M_PI + ((yj > y1) ? 180.0f : 0.0f);
    *theta = atanf(-zj / (machine.y1 - yj)) + ((yj > machine.y1) ? M_PI : 0.0f);

    return 0;
}

// inverse kinematics: (x0, y0, z0) -> (position[X_AXIS], position[Y_AXIS], position[Z_AXIS])
// returned status: 0=OK, -1=non-existing position
int delta_calcInverse (coord_data_t *mpos, float *position)
{
    int status;

    position[X_AXIS] = position[Y_AXIS] = position[Z_AXIS] = 0.0f;

    if((status = delta_calcAngleYZ(mpos->x, mpos->y, mpos->z, &position[X_AXIS])) == 0 &&
        (status = delta_calcAngleYZ(mpos->x * COS120 + mpos->y * SIN120, mpos->y * COS120 - mpos->x * SIN120, mpos->z, &position[Y_AXIS])) == 0)  // rotate coords to +120 deg
        status = delta_calcAngleYZ(mpos->x * COS120 - mpos->y * SIN120, mpos->y * COS120 + mpos->x * SIN120, mpos->z, &position[Z_AXIS]);  // rotate coords to -120 deg

    return status;
}

// Returns machine position in mm converted from system position.
static float *transform_to_cartesian (float *target, float *position)
{
    float y1 = -(machine.t + machine.cfg.rf * cosf(position[X_AXIS]));
    float z1 = -machine.cfg.rf * sinf(position[X_AXIS]);

    float y2 = (machine.t + machine.cfg.rf * cosf(position[Y_AXIS])) * SIN30;
    float x2 = y2 * TAN60;
    float z2 = -machine.cfg.rf * sinf(position[Y_AXIS]);

    float y3 = (machine.t + machine.cfg.rf * cosf(position[Z_AXIS])) * SIN30;
    float x3 = -y3 * TAN60;
    float z3 = -machine.cfg.rf * sinf(position[Z_AXIS]);

    float dnm = (y2 - y1) * x3 -(y3 - y1) * x2;

    float w1 = y1 * y1 + z1 * z1;
    float w2 = x2 * x2 + y2 * y2 + z2 * z2;
    float w3 = x3 * x3 + y3 * y3 + z3 * z3;

    // x = (a1*z + b1)/dnm
    float a1 = (z2 - z1) * (y3 - y1) - (z3 - z1) * (y2 - y1);
    float b1 = -((w2 - w1) * (y3 - y1) - (w3 - w1) *(y2 - y1)) / 2.0f;

    // y = (a2*z + b2)/dnm;
    float a2 = -(z2 - z1) * x3 + (z3 - z1) * x2;
    float b2 = ((w2 - w1)  *x3 - (w3 - w1) * x2) / 2.0;

    // a*z^2 + b*z + c = 0
    float a = a1 * a1 + a2 * a2 + dnm * dnm;
    float b = 2.0f * (a1 * b1 + a2 * (b2 - y1 * dnm) - z1 * dnm * dnm);
    float c = (b2 - y1 * dnm) * (b2 - y1 * dnm) + b1 * b1 + dnm * dnm * (z1 * z1 - machine.re_sqr);

    // discriminant
    float d = b * b - 4.0f * a * c;
    if (d < 0.0f)
        target[X_AXIS] = target[Y_AXIS] = target[Z_AXIS] = NAN; // non-existing point
    else {
        target[Z_AXIS] = -0.5f * (b + sqrtf(d)) / a;
        target[X_AXIS] = (a1 * target[Z_AXIS] + b1) / dnm;
        target[Y_AXIS] = (a2 * target[Z_AXIS] + b2) / dnm;
    }

    return target;
}

// Returns machine position in mm converted from system position steps.
static float *delta_convert_array_steps_to_mpos (float *position, int32_t *steps)
{
    float mpos[N_AXIS];

    uint_fast8_t idx = N_AXIS;
    do {
        idx--;
        mpos[idx] = steps[idx] / settings.axis[idx].steps_per_mm;
    } while(idx);

    return transform_to_cartesian(position, mpos);
}

// Transform absolute position from cartesian coordinate system to delta robot coordinate system
static float *transform_from_cartesian (float *target, float *position)
{
    delta_calcInverse((coord_data_t *)position, target);

    return target;
}

static inline float get_distance (float *p0, float *p1)
{
    uint_fast8_t idx = Z_AXIS + 1;
    float distance = 0.0f;

    do {
        idx--;
        distance += (p0[idx] - p1[idx]) * (p0[idx] - p1[idx]);
    } while(idx);

    return sqrtf(distance);
}

// Delta robots needs long lines divided up
static float *delta_segment_line (float *target, float *position, plan_line_data_t *pl_data, bool init)
{
    static uint_fast16_t iterations;
    static bool segmented;
    static float distance;
    static coord_data_t delta, segment_target, final_target, mpos;

    uint_fast8_t idx = N_AXIS;

    if(init) {

        jog_cancel = false;
        memcpy(final_target.values, target, sizeof(final_target));

        if(delta_calcInverse((coord_data_t *)target, mpos.values) == 0) {

            transform_to_cartesian(segment_target.values, position);

            delta.x = target[X_AXIS] - segment_target.x;
            delta.y = target[Y_AXIS] - segment_target.y;
            delta.z = target[Z_AXIS] - segment_target.z;

            distance = sqrtf(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);

            if((segmented = distance > machine.cfg.sl)) {

                idx = N_AXIS;
                iterations = (uint_fast16_t)ceilf(distance / machine.cfg.sl);

                do {
                    --idx;
                    delta.values[idx] = delta.values[idx] / (float)iterations;
                } while(idx);

            } else {
                iterations = 1;
                memcpy(&segment_target, &final_target, sizeof(coord_data_t));
            }

            distance /= (float)iterations;

        } else { // out of bounds, report? or should never happen?
            iterations = 1;
            memcpy(&final_target, position, sizeof(coord_data_t));
        }

        iterations++; // return at least one iteration

    } else {

        iterations--;

        if(segmented && iterations > 1) {
            do {
                idx--;
                segment_target.values[idx] += delta.values[idx];
            } while(idx);
        } else
            memcpy(&segment_target, &final_target, sizeof(coord_data_t));

        if(delta_calcInverse(&segment_target, mpos.values) != 0) {
            memcpy(&mpos, &machine.last_pos, sizeof(coord_data_t));
            iterations = 0;
        }

        if(!pl_data->condition.rapid_motion && distance != 0.0f) {
            float rate_multiplier = get_distance(mpos.values, machine.last_pos.values) / distance;
            pl_data->feed_rate *= rate_multiplier;
            pl_data->rate_multiplier = 1.0 / rate_multiplier;
        }

        memcpy(&machine.last_pos, &mpos, sizeof(coord_data_t));
    }

    return iterations == 0 || jog_cancel ? NULL : mpos.values;
}

static void get_cuboid_envelope (void)
{
    float maxz = -machine.cfg.e - machine.cfg.f - machine.cfg.re - machine.cfg.rf;
    float minz = -maxz;
    float sr = 1.0f / settings.axis[X_AXIS].steps_per_mm; // Steps/rev -> rad/step, XYZ motors should have the same setting!
    float pos[N_AXIS];
    uint32_t z;
    coord_data_t mpos, home = {
        .x = machine.cfg.home_angle,
        .y = machine.cfg.home_angle,
        .z = machine.cfg.home_angle
    };
    struct {
        float pos[N_AXIS];
    } r[8];

    transform_to_cartesian(mpos.values, home.values);
    machine.home_z = mpos.z;

    // find extents
    for(z = 0; z < settings.axis[X_AXIS].steps_per_mm * 2.0f * M_PI ; ++z) {
        pos[0] = pos[1] = pos[2] = sr * (float)z;
        transform_to_cartesian(mpos.values, pos);
        if(!isnanf(mpos.x)) {
            if(minz > mpos.z)
                minz = mpos.z;
            if(maxz < mpos.z)
                maxz = mpos.z;
        }
    }

    maxz = machine.home_z;
    if(machine.cfg.max_angle != 0.0f) {
        home.x = home.y = home.z = machine.cfg.max_angle;
        transform_to_cartesian(mpos.values, home.values);
        minz = mpos.z;
    }

    float btf = -machine.cfg.b;
    if(minz < btf)
        minz = btf;
    if(maxz < btf)
        maxz = btf;

    bool ok;
    float middlez = (maxz + minz) * 0.5f;
    float original_dist = (maxz - middlez);
    float dist = original_dist * 0.5f;
    float sum = 0.0f;

    machine.min_angle[A_MOTOR] = machine.min_angle[B_MOTOR] = machine.min_angle[C_MOTOR] = 2.0 * M_PI;
    machine.max_angle[A_MOTOR] = machine.max_angle[B_MOTOR] = machine.max_angle[C_MOTOR] = -2.0 * M_PI;

    do {
        sum += dist;
        mpos.x = sum;
        mpos.y = sum;
        mpos.z = middlez + sum;
        if((ok = delta_calcInverse(&mpos, r[0].pos)) == 0) {
            mpos.y = -sum;
            ok = delta_calcInverse(&mpos, r[1].pos) == 0;
        }

        if(ok) {
            mpos.x = -sum;
            ok = delta_calcInverse(&mpos, r[2].pos) == 0;
        }

        if(ok) {
            mpos.y = sum;
            ok = delta_calcInverse(&mpos, r[3].pos) == 0;
        }

        if(ok) {
            mpos.x = sum;
            mpos.z = middlez - sum;
            ok = delta_calcInverse(&mpos, r[4].pos) == 0;
        }

        if(ok) {
            mpos.y = -sum;
            ok = delta_calcInverse(&mpos, r[5].pos) == 0;
        }

        if(ok) {
            mpos.x = -sum;
            ok = delta_calcInverse(&mpos, r[6].pos) == 0;
        }

        if(ok) {
            mpos.y = sum;
            ok = delta_calcInverse(&mpos, r[7].pos) == 0;
        }

        if(!ok) {
            sum -= dist;
            dist *= 0.5f;
        } else {
            uint32_t i;
            for(i = 0; i < 8; ++i) {
                if(machine.min_angle[A_MOTOR] > r[i].pos[A_MOTOR])
                    machine.min_angle[A_MOTOR] = r[i].pos[A_MOTOR];
                if(machine.max_angle[A_MOTOR] < r[i].pos[A_MOTOR])
                    machine.max_angle[A_MOTOR] = r[i].pos[A_MOTOR];
                if(machine.min_angle[B_MOTOR] > r[i].pos[B_MOTOR])
                    machine.min_angle[B_MOTOR] = r[i].pos[B_MOTOR];
                if(machine.max_angle[B_MOTOR] < r[i].pos[B_MOTOR])
                    machine.max_angle[B_MOTOR] = r[i].pos[B_MOTOR];
                if(machine.min_angle[C_MOTOR] > r[i].pos[C_MOTOR])
                    machine.min_angle[C_MOTOR] = r[i].pos[C_MOTOR];
                if(machine.max_angle[C_MOTOR] < r[i].pos[C_MOTOR])
                    machine.max_angle[C_MOTOR] = r[i].pos[C_MOTOR];
            }
        }
    } while(original_dist > sum && dist > 0.1f);

    sys.work_envelope.min.x = -sum;
    sys.work_envelope.min.y = -sum;
    sys.work_envelope.min.z = middlez - sum;
    sys.work_envelope.max.x = sum;
    sys.work_envelope.max.y = sum;
    sys.work_envelope.max.z = middlez + sum - settings.homing.pulloff;

    home.x = home.y = 0;
    home.z = sys.work_envelope.max.z;
    if(delta_calcInverse(&home, pos) == 0)
        machine.cfg.home_angle = pos[A_MOTOR];

    machine.min_angle[A_MOTOR] = machine.min_angle[B_MOTOR] = machine.min_angle[C_MOTOR] = machine.cfg.home_angle;
    if(machine.cfg.max_angle != 0.0f && machine.cfg.max_angle < machine.max_angle[A_MOTOR])
        machine.max_angle[A_MOTOR] = machine.max_angle[B_MOTOR] = machine.max_angle[C_MOTOR] = machine.cfg.max_angle;

    // resolution
    pos[0] = pos[1] = pos[2] = 0.0f;
    transform_to_cartesian(r[0].pos, pos);
    pos[0] = sr;
    transform_to_cartesian(r[1].pos, pos);

    float x = r[0].pos[0] - r[1].pos[0];
    float y = r[0].pos[1] - r[1].pos[1];
    machine.resolution = sqrtf(x * x + y * y); // use as segment length (/ 2)?
    machine.envelope_midz  = middlez;
}

static float *get_homing_target (float *target, float *position)
{
    uint_fast8_t idx = Z_AXIS + 1;

    do {
        idx--;
        if(homing_mode == HomingMode_Pulloff)
            target[idx] = machine.home_pulloff = position[X_AXIS] / machine.resolution / settings.axis[X_AXIS].steps_per_mm;
        else if(homing_mode == HomingMode_Locate)
            target[idx] = position[X_AXIS] / machine.resolution / settings.axis[X_AXIS].steps_per_mm;
        else
            target[idx] = bit_istrue(settings.homing.dir_mask.value, bit(idx)) ? -M_PI : M_PI; // 0.5 revolution
    } while(idx);

    return target;
}

static uint_fast8_t delta_limits_get_axis_mask (uint_fast8_t idx)
{
    return bit(idx);
}

static void delta_limits_set_target_pos (uint_fast8_t idx)
{
    sys.position[idx] = 0;
}

static float homing_cycle_get_feedrate (axes_signals_t cycle, float feedrate, homing_mode_t mode)
{
    homing_mode = mode;
    kinematics.transform_from_cartesian = get_homing_target;

    // TODO: transform feedrate to rad/min

    return feedrate;
}

// Set machine positions for homed limit switches. Don't update non-homed axes.
// NOTE: settings.max_travel[] is stored as a negative value.
static void delta_limits_set_machine_positions (axes_signals_t cycle)
{
//    uint_fast8_t idx = N_AXIS;

    sys.position[X_AXIS] = sys.position[Y_AXIS] = sys.position[Z_AXIS] = (machine.cfg.home_angle + machine.home_pulloff) * settings.axis[X_AXIS].steps_per_mm;

/*
    if(settings.homing.flags.force_set_origin || true) {
        do {
            if(cycle.mask & bit(--idx)) {
                sys.position[idx] = 0;
                sys.home_position[idx] = 0.0f;
            }
        } while(idx);
    } else do {
        if(cycle.mask & bit(--idx)) {
            sys.home_position[idx] = bit_istrue(settings.homing.dir_mask.value, bit(idx))
                                      ? settings.axis[idx].max_travel + pulloff
                                      : - pulloff;
            sys.position[idx] = lroundf(sys.home_position[idx] * settings.axis[idx].steps_per_mm);
        }
    } while(idx);

    */
}

static void delta_go_home (sys_state_t state)
{
    plan_line_data_t plan_data;

    plan_data_init(&plan_data);
    plan_data.condition.rapid_motion = On;

    if(mc_line(sys.home_position, &plan_data)) {
        protocol_buffer_synchronize();
        sync_position();
    }
}

static void delta_homing_complete (bool success)
{
    coord_data_t pos, home = {
        .x = machine.cfg.home_angle,
        .y = machine.cfg.home_angle,
        .z = machine.cfg.home_angle
    };

    kinematics.transform_from_cartesian = transform_from_cartesian;

    if(success) {

        get_cuboid_envelope();
        transform_to_cartesian(pos.values, home.values);
//        sys.work_envelope.max[Z_AXIS] = pos.z;

        sys.home_position[0] = .0f;
        sys.home_position[1] = .0f;
        sys.home_position[2] = machine.cfg.flags.home_to_cuboid_top
                                ? sys.work_envelope.max.z
                                : machine.home_z - settings.homing.pulloff;

        if(machine.cfg.flags.home_to_cuboid_top)
            protocol_enqueue_rt_command(delta_go_home);
    }

    if(on_homing_completed)
        on_homing_completed(success);
}

static void cancel_jog (sys_state_t state)
{
    jog_cancel = true;
}

// Checks and reports if target array exceeds machine travel limits. Returns false if check failed.
static bool is_target_inside_cuboid (float *target, bool is_cartesian)
{
    bool failed = false;
    uint_fast8_t idx = N_AXIS;

    if(is_cartesian && sys.homed.mask) do {
        idx--;
        if(bit_istrue(sys.homed.mask, bit(idx)) && settings.axis[idx].max_travel < -0.0f)
            failed = target[idx] < sys.work_envelope.min.values[idx] || target[idx] > sys.work_envelope.max.values[idx];
    } while(!failed && idx);

    return is_cartesian && !failed;
}

static inline bool pos_ok (coord_data_t *pos)
{
    return !(pos->x < machine.min_angle[A_MOTOR] || pos->x > machine.max_angle[A_MOTOR] ||
              pos->y < machine.min_angle[B_MOTOR] || pos->y > machine.max_angle[B_MOTOR] ||
               pos->z < machine.min_angle[C_MOTOR] || pos->z > machine.max_angle[B_MOTOR]);
}

// Checks and reports if target array exceeds machine travel limits. Returns false if check failed.
static bool delta_check_travel_limits (float *target, bool is_cartesian)
{
    bool failed = false;
    uint_fast8_t idx = N_AXIS;
    coord_data_t pos;

    if(!is_cartesian) {
        if(isnanf(transform_to_cartesian(pos.values, target)[A_MOTOR]))
            return false;
    }

    if(is_target_inside_cuboid(is_cartesian ? target : pos.values, true))
        return true;

    if(machine.cfg.flags.limit_to_cuboid)
        return false;

    if(is_cartesian && (delta_calcInverse((coord_data_t *)target, pos.values) || !pos_ok(&pos)))
        return false;

    if(!is_cartesian)
        memcpy(&pos.values, target, sizeof(coord_data_t));

    if(is_cartesian && sys.homed.mask) do {
        idx--;
        if(bit_istrue(sys.homed.mask, bit(idx)) && settings.axis[idx].max_travel < -0.0f) {
            if(idx > Z_AXIS)
                failed = target[idx] < sys.work_envelope.min.values[idx] || target[idx] > sys.work_envelope.max.values[idx];
            else
                failed = pos.values[idx] < machine.min_angle[idx] || pos.values[idx] > machine.max_angle[idx];
        }
    } while(!failed && idx);

    return !failed;
}

static void delta_apply_jog_limits (float *target, float *position)
{
    if(sys.homed.mask == 0)
        return;

    if(machine.cfg.flags.limit_to_cuboid)
        apply_jog_limits(target, position);

    else if(!is_target_inside_cuboid(target, true)) {

        coord_data_t pos;

//        if(delta_calcInverse((coord_data_t *)target, pos.values)) {

            bool ok = false;
            uint_fast8_t idx = N_AXIS;
            coord_data_t delta;

            do {
                idx--;
                delta.values[idx] = target[idx] - position[idx];
            } while(idx);

            float dist, length = sqrtf(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);

            delta.x /= length;
            delta.y /= length;
            delta.z /= length;
            length = length * 0.5f;
            dist = 0.0;

            do {

                if(ok)
                    length += dist;
                else
                    length -= dist;

                target[X_AXIS] = position[X_AXIS] + delta.x * length;
                target[Y_AXIS] = position[Y_AXIS] + delta.y * length;
                target[Z_AXIS] = position[Z_AXIS] + delta.z * length;

                ok = delta_calcInverse((coord_data_t *)target, pos.values) == 0 && pos_ok(&pos);

                if(dist > machine.resolution)
                    dist *= 0.5f;
                else if(!ok)
                    dist += 0.1f;
                else
                    break;

            } while(!ok || dist > machine.resolution);

            if(!ok)
                memcpy(target, position, sizeof(coord_data_t));

            hal.stream.write("");
//        }
    }
}

static const char *delta_set_axis_setting_unit (setting_id_t setting_id, uint_fast8_t axis_idx)
{
    const char *unit = NULL;

    if(axis_idx <= Z_AXIS) switch(setting_id) {

        case Setting_AxisStepsPerMM:
            unit = "step/rad";
            break;

        case Setting_AxisMaxRate:
            unit = "rad/min";
            break;

        case Setting_AxisAcceleration:
            unit = "rad/sec^2";
            break;

        case Setting_AxisMaxTravel:
        case Setting_AxisBacklash:
            unit = "rad";
            break;

        default:
            break;
    }

    return unit == NULL && on_set_axis_setting_unit != NULL ? on_set_axis_setting_unit(setting_id, axis_idx) : unit;
}

static const setting_group_detail_t kinematics_groups [] = {
    { Group_Root, Group_Kinematics, "Delta robot"}
};

static const setting_detail_t kinematics_settings[] = {
    { Setting_Kinematics0, Group_Kinematics, "Segment length", "mm", Format_Decimal, "0.00", NULL, NULL, Setting_NonCore, &delta_settings.sl, NULL, NULL },
    { Setting_Kinematics1, Group_Kinematics, "Forearm length", "mm", Format_Decimal, "###0.000", NULL, NULL, Setting_NonCore, &delta_settings.re, NULL, NULL },
    { Setting_Kinematics2, Group_Kinematics, "Bicep length", "mm", Format_Decimal, "##0.000", NULL, NULL, Setting_NonCore, &delta_settings.rf, NULL, NULL },
    { Setting_Kinematics3, Group_Kinematics, "Base radius", "mm", Format_Decimal, "##0.000", NULL, NULL, Setting_NonCore, &delta_settings.f, NULL, NULL },
    { Setting_Kinematics4, Group_Kinematics, "End effector radius", "mm", Format_Decimal, "##0.000", NULL, NULL, Setting_NonCore, &delta_settings.e, NULL, NULL },
    { Setting_Kinematics5, Group_Kinematics, "Base to floor", "mm", Format_Decimal, "###0.000", NULL, NULL, Setting_NonCore, &delta_settings.b, NULL, NULL },
    { Setting_Kinematics6, Group_Kinematics, "Home position", "deg", Format_Decimal, "-#0.000", "-90", "90", Setting_NonCore, &delta_settings.home_angle, NULL, NULL },
    { Setting_Kinematics7, Group_Kinematics, "Flags...", NULL, Format_Bitfield, "Home to cuboid top,Limit to cuboid", NULL, NULL, Setting_NonCore, &delta_settings.flags, NULL, NULL },
    { Setting_Kinematics8, Group_Kinematics, "Max angle", "deg", Format_Decimal, "-#0.000", "-90", "135", Setting_NonCore, &delta_settings.max_angle, NULL, NULL },
};

#ifndef NO_SETTINGS_DESCRIPTIONS

static const setting_descr_t kinematics_settings_descr[] = {
    { Setting_Kinematics0, "Size that moves are broken up into to compensate for the non-linear movement of a delta robot." },
    { Setting_Kinematics1, "Forearm length is the length of the connecting rod in millimeters." },
    { Setting_Kinematics2, "Bicep length is the length of the driven arm in millimeters." },
    { Setting_Kinematics3, "Base radius is the radius in millimeters from the center of the drive platform to the bicep drive axle." },
    { Setting_Kinematics4, "End effector radius is the distance from the center of the end effector to the connection points of the main connecting rods." },
    { Setting_Kinematics6, "Home position in degrees from biceps paralell to the floor." },
    { Setting_Kinematics8, "Max travel in degrees from biceps paralell to the floor. Set to 0 to use calculated max." }
};

#endif

static void delta_settings_changed (settings_t *settings, settings_changed_flags_t changed)
{
    memcpy(&machine.cfg, &delta_settings, sizeof(delta_settings_t));

    machine.cfg.home_angle = delta_settings.home_angle * RADDEG;
    machine.cfg.max_angle = delta_settings.max_angle * RADDEG;

    machine.t = (machine.cfg.f - machine.cfg.e) * TAN30_2;
    machine.y1 = -0.5f * TAN30 * machine.cfg.f; // f/2 * tg 30
    machine.y1_sqr = machine.y1 * machine.y1;
    machine.fe = 0.5f * TAN30 * machine.cfg.e;    // shift center to edge

    machine.rf_sqr = machine.cfg.rf * machine.cfg.rf;
    machine.re_sqr = machine.cfg.re * machine.cfg.re;

    get_cuboid_envelope();
}

static void delta_settings_save (void)
{
    hal.nvs.memcpy_to_nvs(nvs_address, (uint8_t *)&delta_settings, sizeof(delta_settings_t), true);
}

static void delta_settings_restore (void)
{
    delta_settings.b = 400.0f;
    delta_settings.e = 24.0f;
    delta_settings.f = 75.0f;
    delta_settings.re = 300.0f;
    delta_settings.rf = 100.0f;
    delta_settings.sl = .2f;
    delta_settings.home_angle = 0.0f;
    delta_settings.flags.mask = 0;

    hal.nvs.memcpy_to_nvs(nvs_address, (uint8_t *)&delta_settings, sizeof(delta_settings_t), true);
}

static void delta_settings_load (void)
{
    if(hal.nvs.memcpy_from_nvs((uint8_t *)&delta_settings, nvs_address, sizeof(delta_settings_t), true) != NVS_TransferResult_OK)
        delta_settings_restore();

    if(on_set_axis_setting_unit == NULL) {
        on_set_axis_setting_unit = grbl.on_set_axis_setting_unit;
        grbl.on_set_axis_setting_unit = delta_set_axis_setting_unit;
    }
}

static setting_details_t setting_details = {
    .groups = kinematics_groups,
    .n_groups = sizeof(kinematics_groups) / sizeof(setting_group_detail_t),
    .settings = kinematics_settings,
    .n_settings = sizeof(kinematics_settings) / sizeof(setting_detail_t),
#ifndef NO_SETTINGS_DESCRIPTIONS
    .descriptions = kinematics_settings_descr,
    .n_descriptions = sizeof(kinematics_settings_descr) / sizeof(setting_descr_t),
#endif
    .load = delta_settings_load,
    .restore = delta_settings_restore,
    .save = delta_settings_save,
    .on_changed = delta_settings_changed
};

static void report_options (bool newopt)
{
    on_report_options(newopt);

    if(!newopt)
        hal.stream.write("[KINEMATICS:Delta v0.02]" ASCII_EOL);
}

static status_code_t delta_info (sys_state_t state, char *args)
{
    hal.stream.write("Delta robot:" ASCII_EOL);
    hal.stream.write(" Rectangular cuboid envelope:" ASCII_EOL "  X: ");
    hal.stream.write(ftoa(sys.work_envelope.min.x, 3));
    hal.stream.write(" to ");
    hal.stream.write(ftoa(sys.work_envelope.max.x, 3));
    hal.stream.write(" mm" ASCII_EOL "  Y: ");
    hal.stream.write(ftoa(sys.work_envelope.min.y, 3));
    hal.stream.write(" to ");
    hal.stream.write(ftoa(sys.work_envelope.max.y, 3));
    hal.stream.write(" mm" ASCII_EOL "  Z: ");
    hal.stream.write(ftoa(sys.work_envelope.min.z, 3));
    hal.stream.write(" to ");
    hal.stream.write(ftoa(sys.work_envelope.max.z, 3));
    hal.stream.write(" mm" ASCII_EOL "  Mid Z:" ASCII_EOL "  ");
    hal.stream.write(ftoa(machine.envelope_midz, 3));
    hal.stream.write(" mm" ASCII_EOL " Resolution:" ASCII_EOL "  ");
    hal.stream.write(ftoa(machine.resolution, 3));
    hal.stream.write(" mm" ASCII_EOL);

    return Status_OK;
}

static const char *delta_setting_get_description (setting_id_t id)
{
    return id >= Setting_AxisStepsPerMM && id <= (Setting_AxisStepsPerMM + Z_AXIS)
            ? "Travel resolution in steps per radian (57.296 degrees)."
            : (on_setting_get_description ? on_setting_get_description(id) : NULL);
}

static void delta_command_help (void)
{
    hal.stream.write("$DELTA - show info about delta robot parameters" ASCII_EOL);

    if(on_report_command_help)
        on_report_command_help();
}

const sys_command_t delta_command_list[] = {
    {"DELTA", delta_info, { .noargs = On } },
};

static sys_commands_t delta_commands = {
    .n_commands = sizeof(delta_command_list) / sizeof(sys_command_t),
    .commands = delta_command_list
};

sys_commands_t *delta_get_commands()
{
    return &delta_commands;
}

// Initialize API pointers for Delta robot kinematics
void delta_robot_init (void)
{
    if((nvs_address = nvs_alloc(sizeof(delta_settings_t)))) {

        kinematics.limits_set_target_pos = delta_limits_set_target_pos;
        kinematics.limits_get_axis_mask = delta_limits_get_axis_mask;
        kinematics.limits_set_machine_positions = delta_limits_set_machine_positions;
        kinematics.transform_from_cartesian = transform_from_cartesian;
        kinematics.transform_steps_to_cartesian = delta_convert_array_steps_to_mpos;
        kinematics.segment_line = delta_segment_line;
        kinematics.homing_cycle_get_feedrate = homing_cycle_get_feedrate;

        grbl.on_jog_cancel = cancel_jog;
        apply_jog_limits = grbl.apply_jog_limits;
        grbl.apply_jog_limits = delta_apply_jog_limits;
        grbl.check_travel_limits = delta_check_travel_limits;

        on_report_options = grbl.on_report_options;
        grbl.on_report_options = report_options;

        delta_commands.on_get_commands = grbl.on_get_commands;
        grbl.on_get_commands = delta_get_commands;

        on_report_command_help = grbl.on_report_command_help;
        grbl.on_report_command_help = delta_command_help;

        on_homing_completed = grbl.on_homing_completed;
        grbl.on_homing_completed = delta_homing_complete;

        on_setting_get_description = grbl.on_setting_get_description;
        grbl.on_setting_get_description = delta_setting_get_description;

        settings_register(&setting_details);
    }
}

#endif
