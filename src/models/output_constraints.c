
#include <math.h>

#include <models/output_constraints.h>

void OutputConstraints_Model196_Hdf5(double* states)
{
    if (states[1] < 1e-12)
        states[1] = 0;
    if (states[2] < 1e-12)
        states[2] = 0;
    if ((states[3] < 1e-12)||(states[3] > 1e200))
        states[3] = 0;
    if (states[4] < 1e-12)
        states[4] = 0;
    else if (states[4] > 1e200)
        states[4] = fmod(states[4], 1e200);
}

void OutputConstraints_Model254_Hdf5(double* states)
{
    if (states[1] < 1e-12)
        states[1] = 0;
    if (states[2] < 1e-12)
        states[2] = 0;
    if (states[3] < 1e-12)
        states[3] = 0;
    if (states[4] < 1e-12)
        states[4] = 0;
    else if (states[4] > 1e200)
        states[4] = fmod(states[4], 1e200);
    if (states[5] < 1e-12)
        states[5] = 0;
    else if (states[5] > 1e200)
        states[5] = fmod(states[5], 1e200);
    if (states[6] < 1e-12)
        states[6] = 0;
}

void OutputConstraints_Model256_Hdf5(double* states)
{
    if (states[1] < 1e-12)
        states[1] = 0;
    if (states[2] < 1e-12)
        states[2] = 0;
    if (states[3] < 1e-12)
        states[3] = 0;
    if (states[4] < 1e-12)
        states[4] = 0;
    else if (states[4] > 1e200)
        states[4] = fmod(states[4], 1e200);
    if (states[5] < 1e-12)
        states[5] = 0;
    else if (states[5] > 1e200)
        states[5] = fmod(states[5], 1e200);
    if (states[6] < 1e-12)
        states[6] = 0;
    else if (states[6] > 1e200)
        states[6] = fmod(states[6], 1e200);
    if (states[7] < 1e-12)
        states[7] = 0;
}