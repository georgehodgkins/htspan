void print_references(plp_aux_t* bfile, size_t n) {
	cout << "reference sequences [" << bfile->hdr->n_targets << "]:" << endl;
	const int total = bfile->hdr->n_targets;
	for (size_t i = 0; i < total; ++i) {
		if (i == n - 2 && n < total) {
			cout << "..." << endl;
		} else if (i > n - 2 && i < bfile->hdr->n_targets - 1) {
			// print nothing
		} else {
			cout << bfile->hdr->target_name[i] << endl;
		}
	}
	cout << endl;
}

void print_n_alignment(plp_aux_t* bfile, size_t n) {
	bam1_t *b = bam_init1();

	cout << "qname\tseq\tqual\tcigar\taux" << endl;

	for (size_t r = 0; r < n; ++r) {

		if (bam_read1(bfile->hf->fp.bgzf, b) > 0) {
			cout << bam_get_qname(b) << '\t';

			for (size_t i = 0; i < b->core.l_qseq; ++i) {
				cout << nuc_to_char(bam_seqi(bam_get_seq(b), i));
			}
			cout << '\t';

			const uint8_t* qual = bam_get_qual(b);
			for (size_t i = 0; i < b->core.l_qseq; ++i) {
				cout << char((*qual) + 33);
				++qual;
			}
			cout << '\t';

			const uint32_t* cigar = bam_get_cigar(b);
			for (size_t i = 0; i < b->core.n_cigar; ++i) {
				cout << *cigar++;
			}
			cout << '\t';

			cout << bam_get_aux(b);

			cout << endl;
		} else {
			cerr << "Error: Could not read record" << endl;
		}

	}
	cout << endl;

	bam_destroy1(b);
}
