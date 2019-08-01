library(io)
library(precrec)
library(ggplot2)
library(jsonlite)


add_id <- function(d) {
	d$snv <- paste(d$chrom, d$pos, d$ref, d$alt, sep="_")
	d
}

annotate_truth <- function(d, positives) {
	d$truth <- d$snv %in% positives$snv;
	d
}

#group <- "01_01_01";
#group <- "01_01_02";
group <- "01_01_03";

out.fname <- filename("eval", path=c("grid", "eval", group), date=NA);
pdf.fname <- insert(out.fname, ext="pdf");

params <- read_json(file.path("grid", group, "params.json"));

positives <- qread(file.path("grid", group, "orig.snv"), type="tsv");
positives <- add_id(positives);

fixed <- qread(file.path("grid", "obias", group, "ffpe.fixed.snv"), type="tsv", fill=TRUE);
fixed <- annotate_truth(add_id(fixed), positives);

unknown <- qread(file.path("grid", "obias", group, "ffpe.unknown.snv"), type="tsv", fill=TRUE);
unknown <- annotate_truth(add_id(unknown), positives);

variable <- qread(file.path("grid", "obias", group, "ffpe.variable.snv"), type="tsv", fill=TRUE);
variable <- annotate_truth(add_id(variable), positives);

summary(fixed$FOBP)
summary(unknown$FOBP)
summary(variable$BOBP)

p.fixed <- with(fixed, evalmod(scores=-FOBP, labels=truth));
p.unknown <- with(unknown, evalmod(scores=-FOBP, labels=truth));
p.variable <- with(variable, evalmod(scores=BOBP, labels=truth));


curvetypes <- c("ROC", "PRC");

get_auc <- function(m, ct) {
	d <- auc(m);
	d$aucs[d$curvetypes == ct]
}

auc_title <- function(m, ct) {
	x <- get_auc(m, ct);
	sprintf("AUC: %.3f", x)
}

options(plot = list(width=2.5, height=2.5, res=300, units="in"));

for (ct in curvetypes) {
	qdraw(
		autoplot(p.fixed, curvetype=ct) + ggtitle(auc_title(p.fixed, ct))
		,
		file=insert(pdf.fname, c("fixed", tolower(ct)))
	);
	qdraw(
		autoplot(p.unknown, curvetype=ct) + ggtitle(auc_title(p.unknown, ct))
		,
		file=insert(pdf.fname, c("unknown", tolower(ct)))
	);
	qdraw(
		autoplot(p.variable, curvetype=ct) + ggtitle(auc_title(p.variable, ct))
		,
		file=insert(pdf.fname, c("variable", tolower(ct)))
	);
}

res <- list(
	mod = list(
		fixed = p.fixed,
		unknown = p.unknown,
		variable = p.fixed
	),
	params = params
);
qwrite(res, file=insert(out.fname, "res", ext="rds"));

