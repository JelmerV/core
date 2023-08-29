/*
  scara.c - scara kinematics implementation

  Part of grblHAL

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

#if SCARA

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "../hal.h"
#include "../settings.h"
#include "../nvs_buffer.h"
#include "../nuts_bolts.h"
#include "../planner.h"
#include "../kinematics.h"
#include "../report.h"
#include "../system.h"

// some config stuff
#define MAX_SEG_LENGTH_MM 2.0f // segmenting long lines due to non-linear motions [mm]

#define A_MOTOR X_AXIS // Lower motor (l1)
#define B_MOTOR Y_AXIS // Upper motor (l2)

// struct to hold the xy coordinates
typedef struct {
    float x;
    float y;
} xy_t;

// struct to hold the joint angles q
typedef struct {
    float q1;
    float q2;
} q_t;

#define ABSOLUTE_JOINT_ANGLES_BIT   bit(1)
#define ELBOW_UP_CONFIGURATION_BIT  bit(2)
#define COUPLE_ON_HOMING_BIT        bit(3)
// struct to hold the machine parameters
typedef struct {
    uint8_t config;
    float l1;
    float l2;
    float home1_offset;
    float home2_offset;
} scara_settings_t;
static scara_settings_t machine = {0}, scara_settings;

// global variables
static bool jog_cancel = false;
static on_report_options_ptr on_report_options;
static on_realtime_report_ptr on_realtime_report;
static on_homing_completed_ptr on_homing_completed;
static nvs_address_t nvs_address;


// ************************ Kinematics Calculations ****************************//

// forward kinematics: (absolute) joint angles to cartesian XY
static xy_t q_to_xy(float q1, float q2) {
    xy_t xy;
    xy.x = machine.l1*cosf(q1*RADDEG) + machine.l2*cosf(q2*RADDEG);
    xy.y = machine.l1*sinf(q1*RADDEG) + machine.l2*sinf(q2*RADDEG);
    return xy;
}

// backwards kinematics: cartesian XY to joint angles (absolute joint angles)
static q_t xy_to_q(float x, float y) {
    q_t q;
    float r_sq = x*x + y*y;
    if (r_sq > (machine.l1 + machine.l2)*(machine.l1 + machine.l2)) {
        q.q1 = q.q2 = NAN;
    } else {
        float cos_q12 = (r_sq - machine.l1*machine.l1 - machine.l2*machine.l2) / (2.0f * machine.l1 * machine.l2);
        float q12 = acosf(cos_q12); //relative angle between l1 and l2
        float beta = atan2f(machine.l2*sinf(q12), machine.l1+machine.l2*cos_q12); //angle between l1 and r

        if (bit_istrue(machine.config, ABSOLUTE_JOINT_ANGLES_BIT)) {
            q.q1 = atan2f(y, x) + beta;
            q12 = -q12;
        } else {
            q.q1 = atan2f(y, x) - beta;
        }

        if (bit_istrue(machine.config, ELBOW_UP_CONFIGURATION_BIT)) {
            q.q2 = q.q1 + q12;
        } else {
            q.q2 = q12;
        }

        //rad to degrees
        q.q1 *= DEGRAD;
        q.q2 *= DEGRAD;
    }
    return q;
}


// finds the max rectangle that fits in the work envelope
static void get_cuboid_envelope (void)
{
    float r_max = machine.l1 + machine.l2;
    float r_min = fabsf(machine.l1 - machine.l2);

    // max rectangle to fit a semicircle
    float max_x = r_max / SQRT2;
    float max_y = r_max / SQRT2;
    
    // using the lower hemisphere. Make config option?
    sys.work_envelope.min[X_AXIS] = -max_x;
    sys.work_envelope.min[Y_AXIS] = -max_y;
    sys.work_envelope.max[X_AXIS] = -r_min;
    sys.work_envelope.max[Y_AXIS] = max_y;
}

// *********************** required grblHAL Kinematics functions ************************ //

// Returns machine position in mm converted from system joint angles
static float *scara_transform_to_cartesian(float *coords, float *angles)
{
    // higher axes unchanged
    uint_fast8_t idx = N_AXIS-1;
    do {
        coords[idx] = angles[idx];
        idx--;
    } while (idx > Y_AXIS);
    
    // apply forward kinematics
    xy_t xy = q_to_xy(angles[A_MOTOR], angles[B_MOTOR]);

    coords[X_AXIS] = xy.x;
    coords[Y_AXIS] = xy.y;

    // char msgOut[100] = {0};
    // snprintf(msgOut, sizeof(msgOut), "[degree2cartesian] q:%0.05f,%0.05f|xy:%.05f,%.05f\n", angles[A_MOTOR], angles[B_MOTOR], xy.x, xy.y);
    // hal.stream.write(msgOut);

    return coords;
}

// Returns machine position in mm converted from system position steps.
static float *scara_transform_steps_to_cartesian(float *position, int32_t *steps)
{
    float angles[N_AXIS] = {0};

    // higher axis dont have to be modified
    uint_fast8_t idx = N_AXIS;
    do {
        idx--;
        angles[idx] = (float)steps[idx] / settings.axis[idx].steps_per_mm;
    } while (idx);

    return scara_transform_to_cartesian(position, angles);
}

// Returns join angles in rad, converted from machine position in mm
static float *scara_transform_from_cartesian(float *target_q, float *position_xy)
{
    // do not change higher axis
    uint_fast8_t idx = N_AXIS-1;
    do {
        target_q[idx] = position_xy[idx];
        idx--;
    } while (idx > Y_AXIS);

    // apply inverse kinematics
    q_t q = xy_to_q(position_xy[A_MOTOR], position_xy[B_MOTOR]);

    // Check if out of reach
    if (isnan(q.q1) || isnan(q.q2)) {
        // trigger soft limit
        // TODO: add smarter way of handling out of reach commands
        system_raise_alarm(Alarm_SoftLimit);
        return NULL;
    }

    // char msgOut[100] = {0};
    // snprintf(msgOut, sizeof(msgOut), "[cartesian2degree] xy:%.05f,%.05f|q:%.05f,%.05f\n", position_xy[X_AXIS], position_xy[Y_AXIS], q.q1, q.q2);
    // hal.stream.write(msgOut);

    target_q[A_MOTOR] = q.q1 ;
    target_q[B_MOTOR] = q.q2 ;

    return target_q;
}


// segment long lines into smaller segments for non-linear kinematics
// target is cartesian, position transformed (joint angles)
// first call: init = true, position is current motor steps, target is cartesian coordinates
// later calls: init = false, position is null, target = init return value, now return: next segment target in joint steps
static float *scara_segment_line (float *target, float *position, plan_line_data_t *plan_data, bool init)
{
    static bool do_segments;
    static uint_fast16_t iterations;
    static coord_data_t delta, segment_target, current_position, final_target;

    uint_fast8_t idx = N_AXIS;
    char msgOut[200] = {0};

    if (init) {
        // resume motion
        jog_cancel = false;

        // save final target
        memcpy(final_target.values, target, sizeof(coord_data_t));

        // get current position in cartesian coordinates
        scara_transform_to_cartesian(current_position.values, position);

        // calculate total delta
        idx = N_AXIS;
        do {
            idx--;
            delta.values[idx] = target[idx] - current_position.values[idx];
        } while(idx);

        // check if segmentation needed
        float distance = sqrtf(delta.x * delta.x + delta.y * delta.y);
        do_segments = !(plan_data->condition.rapid_motion) && distance > MAX_SEG_LENGTH_MM;

        // calculate amount of segments and delta step size
        if (do_segments) {
            iterations = (uint_fast16_t)ceilf(distance / MAX_SEG_LENGTH_MM);
            idx = N_AXIS;
            do {
                idx--;
                delta.values[idx] = delta.values[idx] / (float)iterations;
            } while(idx);

            // save current position as initial segment target
            memcpy(&segment_target, &current_position, sizeof(coord_data_t));
        } 
        else {
            // no segmentation needed: segment target matches final target
            iterations = 1;
            memcpy(&segment_target, &final_target, sizeof(coord_data_t));
        }

        // ensure at least 1 iteration
        iterations++;

        // print debug info
        // snprintf(msgOut, sizeof(msgOut), "seg_line|itrs=%d,do_segments=%d,dist=%f,delta=%f,%f,%f\n", iterations, do_segments, distance, delta.x, delta.y, delta.z);
        // hal.stream.write(msgOut);
    } 
    else {
        // return next segment
        iterations--;
        if(do_segments && iterations > 1) {
            // increment segment target for all axes
            idx = N_AXIS;
            do {
                idx--;
                segment_target.values[idx] += delta.values[idx];
            } while(idx);
        } else {
            // last segment: segment target matches final target
            memcpy(&segment_target, &final_target, sizeof(coord_data_t));
        }
    }
    
    // convert to joint angles
    scara_transform_from_cartesian(current_position.values, segment_target.values);

    // more debug info
    // snprintf(msgOut, sizeof(msgOut), "seg_line|itrs=%d|target_xy=%0.4f,%0.4f|target_q=%0.6f,%0.6f\n", 
    //     iterations, segment_target.x, segment_target.y, current_position.x, current_position.y);
    // hal.stream.write(msgOut);

    if (iterations == 0 || jog_cancel) {
        return NULL;
    } else {
        return current_position.values;
    }
}


static uint_fast8_t scara_limits_get_axis_mask (uint_fast8_t idx)
{   
    if ((idx == A_MOTOR) || (idx == B_MOTOR)) {
        return (bit(X_AXIS) | bit(Y_AXIS));
    } else {
        return bit(idx);
    }
}

static void scara_limits_set_target_pos(uint_fast8_t idx)
{
    q_t q;
    q.q1 = sys.position[X_AXIS] / settings.axis[A_MOTOR].steps_per_mm;
    q.q2 = sys.position[Y_AXIS] / settings.axis[B_MOTOR].steps_per_mm;
    xy_t xy = q_to_xy(q.q1, q.q2);

    char msgOut[100] = {0};
    snprintf(msgOut, sizeof(msgOut), "lim_set_pos: idx:%d|xy:%.05f,%.05f|q:%.05f,%.05f\n", idx, xy.x, xy.y, q.q1, q.q2);
    hal.stream.write(msgOut);

    switch(idx) {
        case X_AXIS:
            sys.position[A_MOTOR] = 0; //q.q1 * settings.axis[A_MOTOR].steps_per_mm;
            break;
        case Y_AXIS:
            sys.position[B_MOTOR] = 0 ; //q.q2 * settings.axis[B_MOTOR].steps_per_mm;
            break;
        default:
            sys.position[idx] = 0;
            break;
    }
}

// Set machine positions for homed limit switches. Don't update non-homed axes.
// NOTE: settings.max_travel[] is stored as a negative value.
static void scara_limits_set_machine_positions (axes_signals_t cycle)
{
    char msgOut[100] = {0};
    snprintf(msgOut, sizeof(msgOut), "scara_limits_set_machine_positions: %d\n", cycle.mask);
    hal.stream.write(msgOut);

    int32_t pulloff = settings.homing.pulloff;
    uint_fast8_t idx = N_AXIS;
    if(settings.homing.flags.force_set_origin || true) {
        do {
            // forced to set pos to 0
            if(cycle.mask & bit(--idx)) {
                sys.position[idx] = 0;
                sys.home_position[idx] = 0.0f;
            }
        } while(idx);
    } else do {
        if(cycle.mask & bit(--idx)) {
            float pulloff = settings.homing.pulloff;
            if (bit_istrue(settings.homing.dir_mask.value, bit(idx)))
                pulloff = -pulloff;
            switch(idx) {
                case A_MOTOR:
                    sys.position[A_MOTOR] = lroundf((machine.home1_offset + pulloff) * settings.axis[A_MOTOR].steps_per_mm);
                    break;
                case B_MOTOR:
                    sys.position[B_MOTOR] = lroundf((machine.home2_offset + pulloff) * settings.axis[B_MOTOR].steps_per_mm);
                    break;
                default:
                    if bit_istrue(settings.homing.dir_mask.value, bit(idx))
                        sys.position[idx] = lroundf((settings.axis[idx].max_travel + settings.homing.pulloff) * settings.axis[idx].steps_per_mm);
                    else
                        sys.position[idx] = lroundf(-settings.homing.pulloff * settings.axis[idx].steps_per_mm);
                    break;
            }
        }
    } while(idx);
}

// custom position transformation for during the homing cycle
static float *get_homing_target (float *target_q, float *position_xy)
{
    // do not change higher axis
    uint_fast8_t idx = N_AXIS-1;
    do {
        target_q[idx] = position_xy[idx];
        idx--;
    } while (idx > Y_AXIS);

    // apply inverse kinematics
    q_t q = xy_to_q(position_xy[A_MOTOR], position_xy[B_MOTOR]);
    
    char msgOut[100] = {0};
    snprintf(msgOut, sizeof(msgOut), "homing_target: %d|xy:%.05f,%.05f|q:%.05f,%.05f\n", idx, position_xy[X_AXIS], position_xy[Y_AXIS], q.q1, q.q2);
    hal.stream.write(msgOut);

    if (bit_istrue(settings.homing.dir_mask.value, bit(A_MOTOR)) || bit_istrue(settings.homing.dir_mask.value, bit(B_MOTOR))) {
        target_q[A_MOTOR] = (machine.home1_offset + HOMING_AXIS_LOCATE_SCALAR*settings.homing.pulloff) * settings.axis[A_MOTOR].steps_per_mm;
        target_q[B_MOTOR] = (machine.home2_offset + HOMING_AXIS_LOCATE_SCALAR*settings.homing.pulloff) * settings.axis[B_MOTOR].steps_per_mm;
    } else {
        target_q[A_MOTOR] = (machine.home1_offset - HOMING_AXIS_LOCATE_SCALAR*settings.homing.pulloff) * settings.axis[A_MOTOR].steps_per_mm;
        target_q[B_MOTOR] = (machine.home2_offset - HOMING_AXIS_LOCATE_SCALAR*settings.homing.pulloff) * settings.axis[B_MOTOR].steps_per_mm;
    }
    // do {
    //     if((settings.homing.pulloff * HOMING_AXIS_LOCATE_SCALAR - fabsf(position_xy[X_AXIS])) > - 0.1f) {
    //         target_q[idx] = position_xy[X_AXIS] * 360.0f; // settings.axis[X_AXIS].steps_per_mm;
    //     } else {
    //         target_q[idx] = bit_istrue(settings.homing.dir_mask.value, bit(idx)) ? -.5f : .5f;
    //     }
    //     idx--;
    // } while(idx);

    return target_q;
}

static float homing_cycle_get_feedrate (float feedrate, axes_signals_t cycle)
{
    char msgOut[100] = {0};
    snprintf(msgOut, sizeof(msgOut), "homing_cycle_get_feedrate: %d\n", cycle.mask);
    hal.stream.write(msgOut);

    // currently homing, do required modifications
    if (bit_istrue(cycle.mask, bit(A_MOTOR)) || bit_istrue(cycle.mask, bit(B_MOTOR))) {      
        kinematics.transform_from_cartesian = get_homing_target;
        
        if (bit_istrue(machine.config, COUPLE_ON_HOMING_BIT)) {
            // couple A and B motors to move in sync
            // ... TODO
            hal.stream.write("couple motors\n");
        }
    }

    return feedrate;
}

static void scara_homing_complete(bool success)
{
    hal.stream.write("scara_homing_complete\n");
    
    kinematics.transform_from_cartesian = scara_transform_from_cartesian;

    if (bit_istrue(machine.config, COUPLE_ON_HOMING_BIT)) {
        // restore coupled motors
        // ... TODO
    }

    if(success) {
        // update work envelope
        get_cuboid_envelope();
    }

    if(on_homing_completed) {
        on_homing_completed(success);
    }
}


// ************************ Settings ************************ //

static const setting_group_detail_t kinematics_groups [] = {
    { Group_Root, Group_Kinematics, "Scara robot"}
};

static const setting_detail_t kinematics_settings[] = {
    { Setting_Kinematics0, Group_Kinematics, "Kinematics configuration", NULL, Format_Bitfield, "Absolute angles, Elbow Up, Couple motors during homing", NULL, NULL, Setting_NonCore, &scara_settings.config, NULL, NULL },
    { Setting_Kinematics1, Group_Kinematics, "Link 1 length", "mm", Format_Decimal, "###0.00", NULL, NULL, Setting_NonCore, &scara_settings.l1, NULL, NULL },
    { Setting_Kinematics2, Group_Kinematics, "Link 2 length", "mm", Format_Decimal, "###0.00", NULL, NULL, Setting_NonCore, &scara_settings.l2, NULL, NULL },
    { Setting_Kinematics3, Group_Kinematics, "Homing switch 1 position", "degrees", Format_Decimal, "##0.00", NULL, NULL, Setting_NonCore, &scara_settings.home1_offset, NULL, NULL },
    { Setting_Kinematics4, Group_Kinematics, "Homing switch 2 position", "degrees", Format_Decimal, "##0.00", NULL, NULL, Setting_NonCore, &scara_settings.home2_offset, NULL, NULL },
};

static void scara_settings_changed (settings_t *settings, settings_changed_flags_t changed)
{
    float position[N_AXIS] = {0}, cartesian[N_AXIS];

    memcpy(&machine, &scara_settings, sizeof(scara_settings_t));

    get_cuboid_envelope();
}

static void scara_settings_save(void)
{
    hal.nvs.memcpy_to_nvs(nvs_address, (uint8_t *)&scara_settings, sizeof(scara_settings_t), true);
}

static void scara_settings_restore(void)
{
    // setting defaults
    scara_settings.config = 0b111;
    scara_settings.l1 = 500.0f;
    scara_settings.l2 = 450.0f;
    scara_settings.home1_offset = 45.0f;
    scara_settings.home2_offset = -180.0f;

    hal.nvs.memcpy_to_nvs(nvs_address, (uint8_t *)&scara_settings, sizeof(scara_settings_t), true);
}

static void scara_settings_load(void)
{
    if(hal.nvs.memcpy_from_nvs((uint8_t *)&scara_settings, nvs_address, sizeof(scara_settings_t), true) != NVS_TransferResult_OK)
        scara_settings_restore();
}

static setting_details_t setting_details = {
    .groups = kinematics_groups,
    .n_groups = sizeof(kinematics_groups) / sizeof(setting_group_detail_t),
    .settings = kinematics_settings,
    .n_settings = sizeof(kinematics_settings) / sizeof(setting_detail_t),
    .load = scara_settings_load,
    .restore = scara_settings_restore,
    .save = scara_settings_save,
    .on_changed = scara_settings_changed
};

// ********************* other functions ********************** //

static void cancel_jog (sys_state_t state)
{
    jog_cancel = true;
}

static void report_options (bool newopt)
{
    on_report_options(newopt);  // call original report before adding new info
    if(!newopt) {
        hal.stream.write("[KINEMATICS:Scara v0.01]" ASCII_EOL);
    }
}

static void report_angles (stream_write_ptr stream_write, report_tracking_flags_t report)
{
    stream_write("|Qj:");
    stream_write(ftoa(sys.position[A_MOTOR]/settings.axis[A_MOTOR].steps_per_mm, 2));
    stream_write(",");
    stream_write(ftoa(sys.position[B_MOTOR]/settings.axis[B_MOTOR].steps_per_mm, 2));
    
    if (on_realtime_report){
        on_realtime_report(stream_write, report);
    }
}

void scara_init_system_position (void)
{
    // set initial position
    sys.position[A_MOTOR] = 0.0;
    if (bit_istrue(machine.config, ELBOW_UP_CONFIGURATION_BIT)) {
        sys.position[B_MOTOR] = lroundf(-90.0f * settings.axis[B_MOTOR].steps_per_mm); //90 degrees down
    } else {
        sys.position[B_MOTOR] = lroundf(90.0f * settings.axis[B_MOTOR].steps_per_mm); //90 degrees up
    }
    hal.stream.write("scara_init_system_position\n");
}


// Initialize API pointers for scara kinematics
void scara_init(void)
{
    if((nvs_address = nvs_alloc(sizeof(scara_settings_t)))) {
        // specify custom kinematics functions
        kinematics.transform_steps_to_cartesian = scara_transform_steps_to_cartesian;
        kinematics.transform_from_cartesian = scara_transform_from_cartesian;
        kinematics.segment_line = scara_segment_line;
        kinematics.init_system_position = scara_init_system_position;
        
        // specify custom homing functions
        kinematics.limits_get_axis_mask = scara_limits_get_axis_mask;
        kinematics.limits_set_target_pos = scara_limits_set_target_pos;
        kinematics.limits_set_machine_positions = scara_limits_set_machine_positions;
        kinematics.homing_cycle_get_feedrate = homing_cycle_get_feedrate;
        
        on_homing_completed = grbl.on_homing_completed;
        grbl.on_homing_completed = scara_homing_complete;

        // jog cancel interrupts line segmentation
        grbl.on_jog_cancel = cancel_jog;

        // add additional report info
        on_report_options = grbl.on_report_options;
        grbl.on_report_options = report_options;

        on_realtime_report = grbl.on_realtime_report;
        grbl.on_realtime_report = report_angles;

        // add scara settings
        settings_register(&setting_details);
    }
}

#endif
