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

#dset <- "grid1";
#groups <- c("01_01_01", "01_01_02", "01_01_03", "01_01_04", "01_01_05");
#groups <- c("05_01_01", "05_01_02", "05_01_03", "05_01_04", "05_01_05");

#dset <- "grid2";
#groups <- c("05_01_01", "05_01_02", "05_01_03", "05_01_04");

#dset <- "grid3";
#groups <- c("01_01_01", "01_01_02", "01_01_03", "01_01_04");
#groups <- c("05_01_01", "05_01_02", "05_01_03", "05_01_04");

dset <- "grid4";
#groups <- c("01_01_01", "01_01_02", "01_01_03", "01_01_04");
groups <- c("05_01_01", "05_01_02", "05_01_03", "05_01_04");

for (group in groups) {

	out.fname <- filename("eval", path=c(dset, "eval", group), date=NA);
	pdf.fname <- insert(out.fname, ext="pdf");

	params <- read_json(file.path(dset, "input", group, "params.json"));

	positives <- qread(file.path(dset, "input", group, "orig.snv"), type="tsv");
	positives <- add_id(positives);

	vaf <- qread(file.path(dset, "vaf", group, "ffpe.vaf.snv"), type="tsv", fill=TRUE);
	vaf <- annotate_truth(add_id(vaf), positives);
	
	sob <- qread(file.path(dset, "sobdetector", group, "ffpe.sobdetector.snv"), type="tsv", fill=TRUE, na=".");
	sob <- annotate_truth(add_id(sob), positives);

	gatk <- qread(file.path(dset, "gatk", group, "ffpe.gatk.snv"), type="tsv", fill=TRUE);
	gatk <- annotate_truth(add_id(gatk), positives);

	# OBP appears to be NA for non-filtered variants; therefore, set NA to 0
	gatk$OBP[is.na(gatk$OBP)] <- 0;
	# however, doing this obliterates GATK's performance, so we ultimately don't do it.

	# remove variants with multiple alternatives
	#alleles <- levels(vaf$ref);
	#gatk <- gatk[gatk$alt %in% alleles & gatk$ref %in% alleles, ];
	#gatk$ref <- factor(gatk$ref, alleles);
	#gatk$alt <- factor(gatk$alt, alleles);

	fixed <- qread(file.path(dset, "obias", group, "ffpe.fixed.snv"), type="tsv", fill=TRUE);
	fixed <- annotate_truth(add_id(fixed), positives);

	unknown <- qread(file.path(dset, "obias", group, "ffpe.unknown.snv"), type="tsv", fill=TRUE);
	unknown <- annotate_truth(add_id(unknown), positives);

	variable <- qread(file.path(dset, "obias", group, "ffpe.variable.snv"), type="tsv", fill=TRUE);
	variable <- annotate_truth(add_id(variable), positives);

	summary(vaf$VAFF)
	summary(sob$SOB)
	summary(gatk$OBP)
	summary(fixed$FOBP)
	summary(unknown$FOBP)
	summary(variable$BOBP)

	p.vaf <- with(vaf, evalmod(scores=-VAFF, labels=truth));
	p.sob <- with(sob, evalmod(scores=-abs(SOB), labels=truth));
	# higher is better
	p.gatk <- with(gatk, evalmod(scores=OBP, labels=truth));
	# lower is better
	p.fixed <- with(fixed, evalmod(scores=-FOBP, labels=truth));
	# lower is better
	p.unknown <- with(unknown, evalmod(scores=-FOBP, labels=truth));
	# higher is better
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
	options(plot.device=NULL);

	for (ct in curvetypes) {
		qdraw(
			autoplot(p.vaf, curvetype=ct) + ggtitle(auc_title(p.vaf, ct))
			,
			file=insert(pdf.fname, c("vaf", tolower(ct)))
		);

		qdraw(
			autoplot(p.sob, curvetype=ct) + ggtitle(auc_title(p.sob, ct))
			,
			file=insert(pdf.fname, c("sobdetector", tolower(ct)))
		);

		qdraw(
			autoplot(p.gatk, curvetype=ct) + ggtitle(auc_title(p.gatk, ct))
			,
			file=insert(pdf.fname, c("gatk", tolower(ct)))
		);

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
			vaf = p.vaf,
			sob = p.sob,
			gatk = p.gatk,
			fixed = p.fixed,
			unknown = p.unknown,
			variable = p.fixed
		),
		params = params
	);
	qwrite(res, file=insert(out.fname, "res", ext="rds"));

}

