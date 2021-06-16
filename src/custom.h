#ifndef CUSTOM_H
#define CUSTOM_H

// This file explains and defines how the custom averager
// works. At a glance, by defining CUSTOM_AVG=1, igrmonty
// will run in "average" mode and output a histogram that
// corresponds to a variable "QTY" that is defined in the 
// sample_zone_photon routine within utils.c.
#define CUSTOM_AVG (0)

// Setting these defines the energy limits for photons we
// should track. Values should be given in Hz.
#define CA_MIN_FREQ (200.e9)
#define CA_MAX_FREQ (300.e9)

// Setting these defines the histogram qualities.
#define CA_MIN_QTY (0.001)
#define CA_MAX_QTY (2000.)
#define CA_NBINS   (200)

// In addition to setting the above values, you will need
// to set QTY0 in sample_zone_photon in the utils.c file!


#endif // CUSTOM_H
