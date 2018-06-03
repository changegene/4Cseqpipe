args <- commandArgs(TRUE)

#inputs for the script:

#data input
tab <- as.character(args[1])

#the size (number of bins!) of the smoothing window
window  <- as.numeric(args[2])

#plot limits
mincoord <- as.numeric(args[3])
maxcoord <- as.numeric(args[4])
abs_max_y = as.numeric(args[5])

#ids of experiments (should be present in the table!)
foc_ids <- eval(parse(text=paste('list(', args[6], ')')))

do_plot = 0
do_table = 1

#output name
tab_nm = as.character(args[7])

#if this is 1, only non-unique fragends will be used, otherwise only unique are considered
multi = as.numeric(args[8])

quants = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)

#For reference
beta_ids=c(498,499,477,524,525,546,489,490,510,511,530,550,493,513,532,533,553,554)
alpha_ids=c(487,488,508,509,529,549,491,492,512,531,551,552)
satb_ids = c(514,515,534,535,555,556)

bd = read.table(tab, header=T)

#to save time, the moving quantile jumps in steps (but we set the step so that our window will have about 40 steps - ensuring high enough resolution
step = round(window/40)

have_ref = 0

rpt_n = length(foc_ids)
nplts = length(foc_ids)*3+1

combined = c()
combined_ord = c()

mat = data.frame()

for(exp_id in foc_ids) {

#extract only the columns of the relevant experiment
	bd_id = bd[,grep(paste(exp_id,"_", sep=""), names(bd))]

#change column names for simplicity
	names(bd_id) = sub(paste("e", exp_id,"_",sep=""), "e", names(bd_id))

#extract 3 and 5 frag end data, unique, non-blind
	uni_cut_3 = bd_id$e3
	uni_cut_3[bd_id$efe3 > bd_id$efl] = NA

	uni_cut_5 = bd_id$e5
	uni_cut_5[bd_id$efe5 > bd_id$efl] = NA
	if(multi == 1) {
		uni_cut_3[bd_id$e3m == 1] = NA
		uni_cut_5[bd_id$e5m == 1] = NA
	} else {
		uni_cut_3[bd_id$e3m != 1] = NA
		uni_cut_5[bd_id$e5m != 1] = NA
	}

#complete missing data with approx (to ensure quantile normalization is not affected by variable density of blind/non_blind or different cutters)
	uni_cut3_full = approx(bd$start, uni_cut_3, n=length(bd$start))
	uni_cut5_full = approx(bd$start, uni_cut_5, n=length(bd$start))

#currenty our reference for the normalization is the first exp in the list
	if(have_ref == 0) {
		reference = c(uni_cut3_full$y, uni_cut5_full$y)
		have_ref = 1
	}

#copmute coumulative density function
	uni_cut3_ecdf = ecdf(uni_cut3_full$y)
	uni_cut5_ecdf = ecdf(uni_cut5_full$y)
#quant normalize given the reference. Not that we put the NA's back (ecdf(uni_cut_3) and not ecdf(uni_cut3_full), but the ecdf is normalized for fragend density

	uni_cut3_norm = quantile(reference, uni_cut3_ecdf(uni_cut_3))
	uni_cut5_norm = quantile(reference, uni_cut5_ecdf(uni_cut_5))

#combining the 3 and 5 ends

	uni_bli_3 = bd_id$e3
	uni_bli_3[bd_id$efe3 < bd_id$efl] = NA

	uni_bli_5 = bd_id$e5
	uni_bli_5[bd_id$efe5 < bd_id$efl] = NA

	if(multi == 1) {
		uni_bli_3[bd_id$e3m == 1] = NA
		uni_bli_5[bd_id$e5m == 1] = NA
	} else {
		uni_bli_3[bd_id$e3m != 1] = NA
		uni_bli_5[bd_id$e5m != 1] = NA
	}

#complete missing data with approx?
	uni_bli3_full = approx(bd$start, uni_bli_3, n=length(bd$start))
	uni_bli5_full = approx(bd$start, uni_bli_5, n=length(bd$start))
#copmute quants
	uni_bli3_ecdf = ecdf(uni_bli3_full$y)
	uni_bli5_ecdf = ecdf(uni_bli5_full$y)
#quant normalize to a reference
	uni_bli3_norm = quantile(reference, uni_bli3_ecdf(uni_bli_3))
	uni_bli5_norm = quantile(reference, uni_bli5_ecdf(uni_bli_5))

	if(do_table == 1) {
		if(dim(mat)[1] == 0) {
			mat = cbind(bd$start, uni_cut3_norm, uni_cut5_norm, uni_bli3_norm, uni_bli5_norm)
		} else {
			mat = cbind(mat, uni_cut3_norm, uni_cut5_norm, uni_bli3_norm, uni_bli5_norm)
		}
	}
}

if(do_table == 1){
	write.table(mat, tab_nm, sep="\t", quote=F,row.names=F,col.names=F)
}

