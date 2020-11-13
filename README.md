# GFU Online vs. Offline processing

This repository explores whether or not there is a noticeable difference in analysis level quantities when comparing the quicker GFU online processing vs. the more sophisticated settings required for the offline processing

## Time-integrated sensitivity
The relevant script to calculate the 8 year sensitivity for the online and offline samples is `time_integrated.py`. The sensitivities for the offline and online samples are comparable, with offline performing about 5% better than online, as seen below:

![Time-integrated sensitivity comparison](/dump/time_integrated_online_vs_offline_comparison.png)

## Transient sensitivity
For shorter timescales, the more advanced reconstructions make less of a difference, the sensitivity is more reflective of the effective area. The sensitivities, calculated in `transient.py`, for the offline and online sample for an injected spectral index of -2 is displayed below:

![Transient sensitivity comparison](/dump/transient_online_vs_offline_comparison_gamma_2.0.png)