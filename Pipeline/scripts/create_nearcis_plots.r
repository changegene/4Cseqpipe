suppressPackageStartupMessages(library(zoo))
require("zoo", quiet=T)

args <- commandArgs(TRUE)
fn <- args[1]
res = as.numeric(args[2])/1000*3-1
feat_fn <- args[3]
message("feat_fn is ", feat_fn)

stat <- args[4]

adjust_color = 1
if(!is.na(args[5])) {
	adjust_color <- as.numeric(args[5])
}

out_fn <- args[6]

color_def_script <- args[7]

plot_from_shift <- args[8]
plot_to_shift <- args[9]

int_type <- args[10]
int_down <- args[11]
int_up <- args[12]
trunc_p <- args[13]
trunc_n <- args[14]
ymax <- as.numeric(args[15])

spec = colorRampPalette(c("white", "gray", "cyan", "lightblue", "blue", "pink", "red"))
spec_func = colorRamp(c("white", "gray", "cyan", "lightblue", "blue", "pink", "red"))
shades = spec(300)

trend_line_color = "black"
trend_line_width = 3
interval_background_color = "lightgray"
points_color = "gray65"

if(!is.na(color_def_script) & color_def_script != "NA") {
	source(color_def_script);
}

system(sprintf("perl scripts/gen_nearcis_scales.pl %s %s %s %s %s %s %s %s %s", fn, stat, stat, res, int_type, int_down, int_up, trunc_p, trunc_n))

sc = read.table(sprintf("%s.%s.scales", fn, stat))
factor_t  = read.table(sprintf("%s.%s.normfactor_t", fn, stat))[1,1]
factor_m  = read.table(sprintf("%s.%s.normfactor_m", fn, stat))[1,1]

factor_t = factor_t*adjust_color

if (feat_fn != "NA"){
	feat = read.table(feat_fn,col.names=c("from", "to", "name", "color"), header=T)
}

png(out_fn, width=2000, height=430)

par(oma=c(0, 2.00, 0, 0.55))
par(mar=c(0, 2.00, 0.5, 0.55))
par(yaxs="i")
par(xaxs="i")
sc[is.na(sc[,res]), res] = 0
sc[is.na(sc[,res+1]), res+1] = 0
sc[is.na(sc[,res+2]), res+2] = 0

data = read.table(fn)

up_smooth = as.vector(rollmean(zoo(sc[,res+1],sc[,1]),7,na.pad=T,fill=T))
down_smooth = as.vector(rollmean(zoo(sc[,res+2],sc[,1]),7,na.pad=T,fill=T))
#up_smooth = as.vector(rollmean(zoo(sc[,res+1],sc[,1]),7,fill=T))
#down_smooth = as.vector(rollmean(zoo(sc[,res+2],sc[,1]),7,fill=T))
up_smooth[is.na(up_smooth)] = 0
down_smooth[is.na(down_smooth)] = 0

#smooth the center very delicately 
center_smooth = as.vector(rollmean(zoo(sc[,res],sc[,1]),3,na.pad=T,fill=T))
#center_smooth = as.vector(rollmean(zoo(sc[,res],sc[,1]),3,fill=T))
center_smooth[is.na(center_smooth)] = 0
center_smooth[1] = 0
center_smooth[length(center_smooth)] = 0

plot_from_shift = as.numeric(plot_from_shift)
plot_to_shift = as.numeric(plot_to_shift)

xfrom = range(sc[,1])[1] + plot_from_shift
xto = range(sc[,1])[2] - plot_to_shift

layout(matrix(c(1,2,3), 3, 1), heights=c(2.5, 3.6, 4.5))
par(mar=c(0, 2.00, 0, 0.55))
plot(sc[,1], center_smooth, type="l", lwd=2, col="white", axes=F, ylim=c(0, 1), xlim=c(xfrom,xto), ylab="")
if (exists("feat") ){
	text((feat[,1]+feat[,2])/2, rep(0.05, length(feat[,1])), labels= feat[,3], col="black", srt=90, cex=1.5, adj=c(0, 0.5))
}

plot(sc[,1], center_smooth, type="l", lwd=2, col="white", ylim=c(0,ymax), xlim=c(xfrom,xto), yaxt='n', xaxt='n', bty='n')
#axis(2, pos=sc[1,1]-(4300/1300000)*(xto-xfrom)+plot_from_shift, lwd=1, col="black", las=2, cex.axis=1.35)

polygon(c(sc[,1],rev(sc[,1])),c(up_smooth, rev(down_smooth)), col=interval_background_color,border=NA)

for (col in 2:dim(data)[2]) {
	points(data[,1], pmin(1.2,data[,col]*factor_t), col=points_color, pch=19,cex=0.5)
}

lines(sc[,1], center_smooth, type="l", lwd=trend_line_width, col=trend_line_color)

if (exists("feat") ){
	points((feat[,1]+feat[,2])/2, rep(ymax, length(feat[,1])), col=feat[,4], cex=1.8, pch=25,bg=feat[,4])
}

par(mar=c(4, 2.00, 0.08, 0.55))
par(mgp=c(2.5, 0.3, 0))

image(as.matrix(log2(0.0005+factor_m/factor_t*sc[(1+plot_from_shift/1000):(dim(sc)[1]-plot_to_shift/1000), 3*(45:2)-4])), col=shades, zlim=c(-11.2,0.2),axes=F, xlab = "Chromosomal Coordinate (Mb)", cex.lab = 2.2)

winsize <-  10^nchar(round(diff(c(xfrom, xto))/2))/10
x.axis <- table(floor(sc[(1+plot_from_shift/1000):(dim(sc)[1]-plot_to_shift/1000), 1]/(winsize)))
x.axis <- cumsum(x.axis) -x.axis +1
x.axis <- x.axis[names(x.axis) !="0"] ## remove zero!
par(cex = 1) 
axis(1, 
     at=x.axis/( length(sc[,1]) - plot_from_shift/1000 - plot_to_shift/1000 ), ## make sequence of 0 till 1 by number of steps of rows
     labels=as.numeric(names(x.axis))*winsize/1e6, padj=0.5)
graphics.off()




