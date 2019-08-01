library(io)
library(ggplot2)
library(jsonlite)
library(reshape2)


groups <- c("01_01_01", "01_01_02", "01_01_03", "01_01_04", "01_01_05");
#groups <- c("05_01_01", "05_01_02", "05_01_03", "05_01_04", "05_01_05");

out.fname <- filename("eval", path=c("grid", "eval-quantify"), date=NA);
pdf.fname <- insert(out.fname, ext="pdf");

truths <- unlist(lapply(groups,
	function(group) {
		params <- read_json(file.path("grid", "input", group, "params.json"));
		params$phi
	}
));

obquants <- unlist(lapply(groups,
	function(group) {
		obquant <- read_json(file.path("grid", "obias", group, "ffpe.obquant.json"));
		obquant$estimate$phi
	}
));

bobquants <- unlist(lapply(groups,
	function(group) {
		bobquant <- read_json(file.path("grid", "obias", group, "ffpe.bobquant.json"));
		with(bobquant, alpha_phi / (alpha_phi + beta_phi))
	}
));

gatks <- unlist(lapply(groups,
	function(group) {
		gatk <- read.table(file.path("grid", "gatk", group, "ffpe.pre_adapter_summary_metrics"), header=TRUE, sep="\t");
		gatk.q <- gatk$TOTAL_QSCORE[which(gatk$ARTIFACT_NAME == "Deamination")];
		10^(-gatk.q/10)
	}
));

d <- data.frame(
	truth = truths,
	htspan_freq = obquants,
	htspan_bayes = obquants,
	gatk = gatks
);

dm <- melt(d, id.vars="truth");

qdraw(
	ggplot(dm, aes(x=truth, y=value, colour=variable, shape=variable, linetype=variable)) +
		geom_abline(intercept=0, slope=1, colour="firebrick", alpha=0.1, size=1.5) +
		geom_abline(intercept=0, slope=0.5, colour="navyblue", alpha=0.1, size=1.5) +
		geom_point(alpha=0.5) + geom_line() + theme_bw() + scale_colour_viridis_d(begin=0.9, end=0.1) +
		xlab("truth") + ylab("estimate") + ggtitle("Damage quantification") + xlim(0, 0.5) + ylim(0, 0.5) +
		theme(legend.title=element_blank())
	,
	width = 4, height = 3,
	file = pdf.fname
)
