# Simple EEG Metrics

<!-- #### A Note: -->
> I created this library to calculate the basic EEG frequency-based metrics.

## Initializing the Library
Load the library by copying or downloading the files to MATLAB and add path to the library using the `addpath(...)` function.

## List of Functions:
<table border="0">
<tr>
    <td> 1. </td>
    <td> <a href = "calculate_eeg_bandpow.m">calculate_eeg_bandpow</a> </td>
    <td> A wrapper over the spectrogram function of the Signal Processing Toolbox, it outputs the frequency-spectra and the band-powers. </td>
</tr><tr>
    <td> 2. </td>
    <td> <a href = "calculate_eeg_laa.m">calculate_eeg_laa</a> </td>
    <td> Calculates the Lateral alpha asymmetry using a list of channel locations. </td>
</tr><tr>
    <td> 3. </td>
    <td> <a href = "calculate_eeg_tbr.m">calculate_eeg_tbr</a> </td>
    <td> Calculates the Theta-to-Beta bandpower ratio. </td>
</tr><tr>
    <td> 4. </td>
    <td> <a href = "calculate_eeg_bar.m">calculate_eeg_bar</a> </td>
    <td> Calculates the Beta-to-Alpha bandpower ratio. </td>
</tr><tr>
    <td> 5. </td>
    <td> <a href = "calculate_eeg_metrics.m">calculate_eeg_metrics</a> </td>
    <td> A wrapper over calculate_eeg_laa, calculate_eeg_tbr, and calculate_eeg_bar. </td>
</tr>
</table>
<!-- <p style="font-size:5px"><i>more to follow...</i></p> -->

>Functional dependencies: No external dependencies exist.

# Licensing

The repository is licensed under the terms of [MIT License](LICENSE).
