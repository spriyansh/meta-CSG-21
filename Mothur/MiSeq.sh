# Run 1

make.file(inputdir=., type=fastq, prefix=stability)

make.contigs(file=stability.files)

summary.seqs(fasta=stability.trim.contigs.fasta)

##########################################

# Run 2, choose 267

screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, summary=stability.trim.contigs.summary, maxhomop=8$

summary.seqs(fasta=stability.trim.contigs.good.fasta)

unique.seqs(fasta=stability.trim.contigs.good.fasta)

count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)

summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)

align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.nr_v138_1.align)

summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)


#####################################################

# Run 3, 13144  23961

screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.$

summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)

filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filte$

chimera.vsearch(vsearch=/home/psrivastava/.conda/envs/thesis/bin/vsearch,fasta=stability.trim.contigs.good.unique.good.filter.uniq$

remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.uniqu$

summary.seqs(fasta=current, count=current)

classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good$

remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.goo$

summary.tax(taxonomy=current, count=current)

get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stab$

seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=stability.trim.contig$

# Remove Mock

remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=s$

####################################################

cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=stability.trim.co$

################################################

# OTU

make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.list, count=stability.tri$

classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.list, count=stability.tr$

get.oturep(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist, fasta=stability.trim.contig$

# Complete
