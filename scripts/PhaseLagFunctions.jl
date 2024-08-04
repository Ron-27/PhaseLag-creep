"""
    data_read(filepath)

Creates a DataFrame with the experiment data, based on the structure
of the CSV file.

###### Arguments
- `filepath`: the path to the data file

###### Returns
- `df`: DataFrame with the raw data.

"""
function data_read(filepath)
	df = CSV.read(filepath, DataFrame, header=2)
	rename!(df, "s" => :time_s)
	rename!(df, "mN" => :load_mN)
	rename!(df, "nm" => :eps_nm)
	return df
end


"""

	time_fix(t, round_digit)

Rounding the time signatures according to a number of digits, and calculates the smapling frequency. This helps the FFT and phase calculations.

###### Arguments
- `t`: time vector

- `round_digit`: the number of decimal point digits to round to

###### Returns
- `fs_rounded`: sampling frequency.
- `time_new`: new time vector.
"""
function time_fix(t, round_digit)
	differences = t[2:end] - t[1:end-1]
	Ts_rounded = round(mean(differences), digits=round_digit)
	fs_rounded = 1/Ts_rounded
	start_value = t[1]  # Replace with start value
	step = Ts_rounded         # Replace with step
	len = length(t)      # Replace with length

	# Create the array
	time_new = collect(range(start_value, step=step, length=len))
	
	return fs_rounded, time_new
end


"""

	df_good(df, digit)

Creates a new DataFrame with fixed time signatures and only load and strain data.

###### Arguments
- `df`: DataFrame of original dataset

- `digit`: the number of decimal point digits to round to

###### Returns
- `df_g`: new DataFrame
"""
function df_good(df, digit)
	fs, time = time_fix(df.time_s, digit)
	df_g = df[:, [:load_mN, :eps_nm]]
	df_g[!, :time_s] = time
	df_g=df_g[:,[:time_s, :load_mN, :eps_nm]]
	return df_g
end


"""

	forcing_freq(time, signal)

Calculates the forcing frequency of the experiment from the FFT spectrum

###### Arguments
- `time`: time vector

- `signal`: measured signal

###### Returns
- `f`: forcing frequency
"""
function forcing_freq(time, signal)
	F = rfft(signal .- mean(signal))
	freqs = rfftfreq(length(time), 1/(time[2] - time[1]))
	idx = argmax(abs.(F))
	f = freqs[idx]
	return f
end


"""

	filter_high(signal, f, fs)

High-pass filter that filters out all the frequencis below half of the forcing frequency. This help to clear the trends in the signal if the trends have much lower frequency behavior than the forcing 

###### Arguments
- `signal`: vector of measured signal

- `f`: forcing frequency

- `fs`: sampling frequency

###### Returns
- `filtered_2`: the filtered signal 
"""
function filter_high(signal, f, fs)
	response = Highpass((f/1.5); fs)
	design = Butterworth(2)
	filtered = filt(digitalfilter(response,design), signal)
	filtered_2 = filtered .- mean(filtered)
	return filtered_2
end


"""

	FFT_phase_rel(sig1, sig2, t, fs)

High-pass filter that filters out all the frequencis below half of the forcing frequency. This help to clear the trends in the signal if the trends have much lower frequency behavior than the forcing 

###### Arguments
- `sig1`: vector of measured signal 1

- `sig2`: vector of measured signal 2

- `t`: time vector (assuming the same for both signals)

- `fs`: sampling frequency (assuming the same for both signals)

###### Returns
- `ϕ`: the relative phase angle between the two signals
"""
function FFT_phase_rel(sig1, sig2, t, fs)
	F1 = rfft(sig1)
	F2 = rfft(sig2)
	freqs = rfftfreq(length(t), fs)
	idx = argmax(abs.(F1))
	idx2 = argmax(abs.(F2))
	ϕ1 = angle(F1[idx])
	ϕ2 = angle(F2[idx2])
	ϕ = abs(mod2pi(ϕ1 - ϕ2 + π) - π)
	return ϕ
end


function syn_filter(signal, f, fs)
	response = Highpass((f/2); fs)
	design = Butterworth(2)
	filtered = filt(digitalfilter(response,design), signal)
	filtered_2 = filtered .- mean(filtered)
	return filtered_2
end


function FFT_phase(sig, t, fs)
	F = rfft(sig)
	freqs = rfftfreq(length(t), fs)
	idx = argmax(abs.(F))
	ϕ = -1*(angle(F[idx]) + π/2) # forcing is sin function minus phase
	#ϕ = -angle(F[idx]) # forcing is cos function minus phase

	return ϕ
end