using SpillPap
using Base.Test

i = rand(1:2000)
w = rand(250:600)
l = rand(2:4)
H = rand(200:1000)
run(`Rscript Rbench.R $i $H $w $l`)


# write your own tests here
data = readcsv("data.csv")
data = log(data)

nocorr = fftSpilloverTableDec12(data[(1+i):(w+i),:], l, H, "Const", [pi + 0.001, pi/5, pi/10, 0], true)
corr = fftSpilloverTableDec12(data[(1+i):(w+i),:], l, H, "Const", [pi + 0.001, pi/5, pi/10, 0], false)
tab = readcsv("table_T.txt")
@test all(abs(tab-[nocorr[1][:,:,1]; nocorr[1][:,:,2]; nocorr[1][:,:,3]]*100).<0.1)
run(`rm table_T.txt`)
tab = readcsv("table_F.txt")
@test all(abs(tab-[corr[1][:,:,1]; corr[1][:,:,2]; corr[1][:,:,3]]*100).<0.1)
run(`rm table_F.txt`)
vec = readcsv("overall_FF.txt")
@test all(abs(collect(vec)-collect(corr[2][1])*100).<0.1)
run(`rm overall_FF.txt`)
vec = readcsv("overall_FT.txt")
@test all(abs(collect(vec)-collect(corr[2][2])*100).<0.1)
run(`rm overall_FT.txt`)
vec = readcsv("overall_TF.txt")
@test all(abs(collect(vec)-collect(nocorr[2][1])*100).<0.1)
run(`rm overall_TF.txt`)
vec = readcsv("overall_TT.txt")
@test all(abs(collect(vec)-collect(nocorr[2][2])*100).<0.1)
run(`rm overall_TT.txt`)

# From within sum has to be the same as overal sums
@test collect(corr[2][1])==collect(sum(corr[3][1],1))
# To within sum has to be the same as overal sums
@test all(abs(collect(corr[2][1])-collect(sum(corr[4][1],1))).<1e-10)
# Reconstruction should work
@test sum(collect(corr[2][2]))==spilloverDiebold12(data[(1+i):(w+i),:], l, H, "Const")