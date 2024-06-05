seed = -1
seqfile =  Matrix.phy
treefile = Pori_Cauchy_SM_SpiPori_573.tre
mcmcfile = mcmc.txt
outfile = out.txt
ndata = 5
seqtype = 2 * 0 : nucleotides; 1: codons; 2: AAs
usedata = 0 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 1 * 1: global clock; 2: independent; and 3: correlated rates
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = 1 1 0.1 * birth, death, sampling
rgene_gamma = 2 20 1 * gammaDir prior for rate for genes
sigma2_gamma = 1 10 1 * gammaDir prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...
print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 100000
sampfreq = 100
nsample = 20000
