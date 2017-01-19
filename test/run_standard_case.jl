data=readcsv("/Users/tomaskrehlik/Documents/PHD/2014_Spillovers_final/data_spills.csv")

using SpillPap
using VARmodels

window = 100;
H = 600;
p = 2;
typ = "Const";
bounds = [π + 0.0001; π/5; π/20; 0];

i = 20

est = varEstimate(data[(1+i):(window+i),:], p, typ)
# Getting the fevd tables for the data and testing whether it works
t = spilloverTable(est, H; fevd_function = genFEVD, nocorr = false)
t2 = spilloverTable(est, H; fevd_function = fftGenFEVD, nocorr = false, bounds = bounds)
t3 = spilloverTable(est, H; fevd_function = genFEVD, nocorr = true)
t4 = spilloverTable(est, H; fevd_function = fftGenFEVD, nocorr = false, bounds = bounds)

# Getting the SpilloverTable and SpilloverTableFrequency objects from the estimated tables
SpilloverTable(t)
SpilloverTable(t3)
SpilloverTableFrequency(t2, bounds)
SpilloverTableFrequency(t4, bounds)

# Checking whether the bootstrap procedures work
reps = 300
boot_sample_length=200
burnout=100
fevd_function = genFEVD
nocorr = false

SpilloverTableBooted(est, H, reps; boot_sample_length=boot_sample_length, burnout=burnout, fevd_function=fevd_function, nocorr=nocorr)
SpilloverTableFrequencyBooted(est, H, reps; bounds = bounds, boot_sample_length=boot_sample_length, burnout=burnout, fevd_function=fftGenFEVD, nocorr=nocorr)
