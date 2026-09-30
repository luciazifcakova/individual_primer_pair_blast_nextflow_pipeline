nextflow.enable.dsl=2
include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_CUSTOMER } from './modules/nf-core/multiqc/main'
include { SEQKIT_STATS as SEQKIT_UNFILTERED } from './modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_FILTERED } from './modules/nf-core/seqkit/stats/main'
include { DEMUX_PRECOMPUTE } from './modules/local/demux_precompute'
include { DEMUX_CHUNK } from './modules/local/demux_chunk'
include { MERGE_CHUNKS } from './modules/local/merge_chunks'
include { MERGE_ANALYSIS_SAMPLES } from './modules/local/merge_analysis_samples'
include { PREPARE_DIRECT_INPUTS } from './modules/local/prepare_direct_inputs'
include { FILTER_READS } from './modules/local/filter_reads'
include { PREPARE_QUERIES } from './modules/local/prepare_queries'
include { SPLIT_BLAST_FASTQ } from './modules/local/split_blast_fastq'
include { BLAST_CHUNK } from './modules/local/blast_chunk'
include { BLAST_CLASSIFY } from './modules/local/blast_classify'
include { BLAST_FASTA } from './modules/local/blast_fasta'
include { KRONA_TSV } from './modules/local/krona_tsv'
include { KRONA_HTML } from './modules/local/krona_html'
include { EXCEL_REPORT } from './modules/local/excel_report'
include { MERGE_SEQKIT_STATS } from './modules/local/merge_seqkit_stats'
include { README_DOCX } from './modules/local/readme_docx'
include { PACKAGE_PROJECT } from './modules/local/package_project'

process NORMALIZE_INPUTS {
 label 'process_low'
 input: path sample_map; path project_map; path assay_map; path helper
 output: path 'sample_map.tsv',emit:sample_map; path 'project_map.tsv',emit:project_map; path 'assay_map.tsv',emit:assay_map
 script:"""
 python3 ${helper} ${sample_map} sample_map.tsv
 python3 ${helper} ${project_map} project_map.tsv
 python3 ${helper} ${assay_map} assay_map.tsv
 """
}

process NORMALIZE_DEMUX_PRIMERS {
 label 'process_low'
 input: path primers; path helper
 output: path 'primers.tsv',emit:primers
 script:"""
 python3 ${helper} ${primers} primers.tsv
 """
}

process PREPARE_ANALYSIS_MAP {
 label 'process_low'
 input: path sample_map; path project_map; path assay_map; path helper
 output: path 'analysis_map.tsv',emit:map
 script:"""
 python3 ${helper} --sample-map ${sample_map} --project-map ${project_map} --assay-db-map ${assay_map} --output analysis_map.tsv --validate-databases
 """
}

process SPLIT_FASTQ {
 label 'process_medium'
 input: path fastqs; path helper
 output: path 'chunks/chunk_*.fastq.gz',emit:chunks
 script:"""
 mkdir -p chunks
 python3 ${helper} --reads ${params.demux_chunk_reads} --outdir chunks ${fastqs.join(' ')}
 """
}

def boolean_param(value, name) {
 def text=value?.toString()?.trim()?.toLowerCase()
 if(value instanceof Boolean) return value
 if(text in ['true','1','yes','y']) return true
 if(text in ['false','0','no','n']) return false
 throw new IllegalArgumentException("--${name} must be true or false; received ${value}")
}

workflow {
 do_demultiplex=boolean_param(params.demultiplex,'demultiplex')
 if((!params.raw_fastq && !params.raw_fastq_dir)||!params.sample_map||!params.project_sample_map||!params.amplicon_primers||!params.assay_blast_db_map||!params.taxonomy_db_dir) error 'Missing required parameters: FASTQ input, sample_map, project_sample_map, amplicon_primers, assay_blast_db_map and taxonomy_db_dir are always required'
 if(do_demultiplex && (!params.all_gene_primers||!params.demux_script)) error '--demultiplex true additionally requires all_gene_primers and demux_script'
 raw=Channel.fromPath(params.raw_fastq?:"${params.raw_fastq_dir}/*.{fastq,fq,fastq.gz,fq.gz}",checkIfExists:true).collect()
 sm=Channel.value(file(params.sample_map,checkIfExists:true));pm=Channel.value(file(params.project_sample_map,checkIfExists:true));am=Channel.value(file(params.assay_blast_db_map,checkIfExists:true))
 normh=Channel.value(file("${projectDir}/bin/normalize_table.py"));maph=Channel.value(file("${projectDir}/bin/prepare_analysis_map.py"));splith=Channel.value(file("${projectDir}/bin/split_fastq.py"));blast_h=Channel.value(file("${projectDir}/bin/blast_postprocess.py"));report_h=Channel.value(file("${projectDir}/bin/report_tables.py"));krona_h=Channel.value(file("${projectDir}/bin/krona_tsv.py"));readme_h=Channel.value(file("${projectDir}/bin/readme_docx.py"))
 amplicon_primers=Channel.value(file(params.amplicon_primers,checkIfExists:true));prepare_h=Channel.value(file("${projectDir}/bin/prepare_queries.py"));manifest_h=Channel.value(file("${projectDir}/bin/primer_manifest.py"));split_fasta_h=Channel.value(file("${projectDir}/bin/split_fasta.py"))
 NORMALIZE_INPUTS(sm,pm,am,normh)
 nsm=NORMALIZE_INPUTS.out.sample_map.first();npm=NORMALIZE_INPUTS.out.project_map.first();nam=NORMALIZE_INPUTS.out.assay_map.first()
 PREPARE_ANALYSIS_MAP(nsm,npm,nam,maph)
 def source_reads
 if(do_demultiplex) {
   pr=Channel.value(file(params.all_gene_primers,checkIfExists:true));ds=Channel.value(file(params.demux_script,checkIfExists:true))
   NORMALIZE_DEMUX_PRIMERS(pr,normh); npr=NORMALIZE_DEMUX_PRIMERS.out.primers.first()
   DEMUX_PRECOMPUTE(nsm,npr,ds); models=DEMUX_PRECOMPUTE.out.models.first()
   SPLIT_FASTQ(raw,splith); chunks=SPLIT_FASTQ.out.chunks.flatten().map{f->tuple([id:f.simpleName],f)}
   DEMUX_CHUNK(chunks,models,ds)
   source_parts=DEMUX_CHUNK.out.demux.flatMap{m,d->d.toFile().listFiles().findAll{x->x.isDirectory()&&new File(x,'reads.fastq.gz').isFile()}.collect{x->tuple(x.name,file(new File(x,'reads.fastq.gz')))}}
   MERGE_CHUNKS(source_parts.groupTuple(by:0))
   source_reads=MERGE_CHUNKS.out.reads
 } else {
   direct_h=Channel.value(file("${projectDir}/bin/prepare_direct_inputs.py"))
   PREPARE_DIRECT_INPUTS(nsm,raw,direct_h)
   source_reads=PREPARE_DIRECT_INPUTS.out.reads.flatten().map{f->tuple(f.name.replaceFirst(/(?i)\.fastq\.gz$/,''),f)}
 }
 rows=PREPARE_ANALYSIS_MAP.out.map.splitCsv(header:true,sep:'\t',strip:true)
 joined=rows.map{r->tuple(r.sample,r.project_id,r.analysis_sample,r.assay,r.blast_db,r.min_quality,r.min_length,r.max_length)}.join(source_reads,by:0)
 grouped=joined.map{s,p,a,assay,db,minq,minl,maxl,fq->tuple(p,a,assay,db,minq,minl,maxl,fq)}.groupTuple(by:[0,1,2,3,4,5,6]).map{p,a,assay,db,minq,minl,maxl,fqs->tuple(p,a,fqs)}
 MERGE_ANALYSIS_SAMPLES(grouped)
 info=rows.map{r->tuple("${r.project_id}\t${r.analysis_sample}",[r.assay,r.blast_db,r.min_quality,r.min_length,r.max_length])}.unique()
 analysis=MERGE_ANALYSIS_SAMPLES.out.reads.map{p,a,fq->tuple("${p}\t${a}",fq)}.join(info,by:0).map{k,fq,inf->
   def q=(inf[2]?.toString()?.trim()) ? inf[2].toFloat() : params.min_quality
   def minl=(inf[3]?.toString()?.trim()) ? inf[3].toInteger() : params.min_length
   def maxl=(inf[4]?.toString()?.trim()) ? inf[4].toInteger() : params.max_length
   tuple([id:k.split('\t')[1],project_id:k.split('\t')[0],assay:inf[0],blast_db:inf[1],chopper_min_quality:q,chopper_min_length:minl,chopper_max_length:maxl,control:(k.split('\t')[1]=~params.control_regex).find()],fq)
 }
 SEQKIT_UNFILTERED(analysis); FILTER_READS(analysis)
 PREPARE_QUERIES(FILTER_READS.out.reads,amplicon_primers,prepare_h,manifest_h)
 SEQKIT_FILTERED(PREPARE_QUERIES.out.reads); FASTQC(PREPARE_QUERIES.out.reads)

 // BLAST is split independently for each analysis sample. Every chunk is a separately cached SLURM task.
 SPLIT_BLAST_FASTQ(PREPARE_QUERIES.out.queries.map{m,q,w->tuple(m,q)},split_fasta_h)
 blast_chunks=SPLIT_BLAST_FASTQ.out.chunks.flatMap{m,fs->
   def xs=(fs instanceof List)?fs:[fs]
   xs.collect{f->tuple(m+[chunk_id:f.simpleName],f)}
 }
 BLAST_CHUNK(blast_chunks)
 blast_grouped=BLAST_CHUNK.out.hits
   .map{m,h,c->tuple(m.project_id,m.id,m.assay,m.blast_db,m.control,h,c)}
   .groupTuple(by:[0,1,2,3,4])
   .map{p,s,assay,db,ctrl,hs,cs->tuple([id:s,project_id:p,assay:assay,blast_db:db,control:ctrl],hs,cs)}
 query_weights=PREPARE_QUERIES.out.queries.map{m,q,w->tuple("${m.project_id}\t${m.id}",w)}
 blast_with_weights=blast_grouped.map{m,h,c->tuple("${m.project_id}\t${m.id}",m,h,c)}
   .join(query_weights,by:0).map{k,m,h,c,w->tuple(m,h,c,w)}
 BLAST_CLASSIFY(blast_with_weights,blast_h)

 blast_for_fasta=BLAST_CLASSIFY.out.results.map{m,rawt,filt,assign,hits,ra,fa->tuple(m,ra,fa)}; BLAST_FASTA(blast_for_fasta)
 raw_tables=BLAST_CLASSIFY.out.results.map{m,rawt,filt,assign,hits,ra,fa->tuple(m,'RAW',rawt)}
 filt_tables=BLAST_CLASSIFY.out.results.map{m,rawt,filt,assign,hits,ra,fa->tuple(m,'FILTERED',filt)}
 KRONA_TSV(raw_tables.mix(filt_tables),krona_h)
 // Customer and control/internal reporting are generated in parallel.
 customer_raw0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->!m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,'RAW','customer',rawt,assign)}.groupTuple(by:[0,1,2])
 customer_fil0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->!m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,'FILTERED','customer',filt,assign)}.groupTuple(by:[0,1,2])
 internal_raw0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,'RAW','internal',rawt,assign)}.groupTuple(by:[0,1,2])
 internal_fil0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,'FILTERED','internal',filt,assign)}.groupTuple(by:[0,1,2])

 customer_qc=BLAST_CLASSIFY.out.qc.filter{m,q->!m.control}.map{m,q->tuple(m.project_id,q)}.groupTuple(by:0)
 internal_qc=BLAST_CLASSIFY.out.qc.filter{m,q->m.control}.map{m,q->tuple(m.project_id,q)}.groupTuple(by:0)
 customer_raw=customer_raw0.join(customer_qc,by:0)
 customer_fil=customer_fil0.join(customer_qc,by:0)
 internal_raw=internal_raw0.join(internal_qc,by:0)
 internal_fil=internal_fil0.join(internal_qc,by:0)
 EXCEL_REPORT(customer_raw.mix(customer_fil).mix(internal_raw).mix(internal_fil),report_h)

 customer_krona=KRONA_TSV.out.tsv.filter{m,mode,f->!m.control}.map{m,mode,f->tuple(m.project_id,mode,'customer',f)}.groupTuple(by:[0,1,2])
 internal_krona=KRONA_TSV.out.tsv.filter{m,mode,f->m.control}.map{m,mode,f->tuple(m.project_id,mode,'internal',f)}.groupTuple(by:[0,1,2])
 KRONA_HTML(customer_krona.mix(internal_krona))

 ust=SEQKIT_UNFILTERED.out.stats.filter{m,f->!m.control}.map{m,f->tuple(m.project_id,'FASTQ_unfiltered_statistics',f)}.groupTuple(by:[0,1])
 fst=SEQKIT_FILTERED.out.stats.filter{m,f->!m.control}.map{m,f->tuple(m.project_id,'FASTQ_statistics',f)}.groupTuple(by:[0,1])
 MERGE_SEQKIT_STATS(ust.mix(fst))

 // Separate MultiQC reports for customer samples and NK/PK controls.
 qcf_customer=FASTQC.out.zip.filter{m,f->!m.control}.map{m,f->tuple(m.project_id,f)}.groupTuple(by:0).map{p,fs->tuple([id:"${p}_customer",project_id:p,kind:'customer'],fs)}
 qcf_internal=FASTQC.out.zip.filter{m,f->m.control}.map{m,f->tuple(m.project_id,f)}.groupTuple(by:0).map{p,fs->tuple([id:"${p}_internal",project_id:p,kind:'internal'],fs)}
 MULTIQC_CUSTOMER(qcf_customer.mix(qcf_internal),Channel.value(file(params.multiqc_config)))

 pids=rows.map{it.project_id}.unique(); README_DOCX(pids,readme_h)

 rawrep=EXCEL_REPORT.out.report.filter{p,mode,kind,x,t->mode=='RAW'&&kind=='customer'}.map{p,m,k,x,t->tuple(p,x,t)}
 filrep=EXCEL_REPORT.out.report.filter{p,mode,kind,x,t->mode=='FILTERED'&&kind=='customer'}.map{p,m,k,x,t->tuple(p,x,t)}
 rawhtml=KRONA_HTML.out.html.filter{p,m,k,h->m=='RAW'&&k=='customer'}.map{p,m,k,h->tuple(p,h)}
 filhtml=KRONA_HTML.out.html.filter{p,m,k,h->m=='FILTERED'&&k=='customer'}.map{p,m,k,h->tuple(p,h)}
 rawk=KRONA_TSV.out.tsv.filter{m,mode,f->!m.control&&mode=='RAW'}.map{m,mode,f->tuple(m.project_id,f)}.groupTuple(by:0)
 filtk=KRONA_TSV.out.tsv.filter{m,mode,f->!m.control&&mode=='FILTERED'}.map{m,mode,f->tuple(m.project_id,f)}.groupTuple(by:0)

 // BLAST zips include one decision per read/OTU, all HSPs, QC and preparation provenance.
 rawb0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->!m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,[rawt,assign,hits])}.groupTuple(by:0).map{p,lists->tuple(p,lists.flatten())}
 filtb0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->!m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,[filt,assign,hits])}.groupTuple(by:0).map{p,lists->tuple(p,lists.flatten())}
 blastqc=BLAST_CLASSIFY.out.qc.filter{m,q->!m.control}.map{m,q->tuple(m.project_id,q)}.groupTuple(by:0)
 preparation_customer=PREPARE_QUERIES.out.audit.filter{m,d->!m.control}.map{m,d->tuple(m.project_id,d)}.groupTuple(by:0)
 preparation_internal=PREPARE_QUERIES.out.audit.filter{m,d->m.control}.map{m,d->tuple(m.project_id,d)}.groupTuple(by:0)
 rawb=rawb0.join(blastqc,by:0).join(preparation_customer,by:0).map{p,files,qcs,dirs->tuple(p,files+qcs+dirs)}
 filtb=filtb0.join(blastqc,by:0).join(preparation_customer,by:0).map{p,files,qcs,dirs->tuple(p,files+qcs+dirs)}
 rawfa=BLAST_FASTA.out.fasta.filter{m,r,f->!m.control}.map{m,r,f->tuple(m.project_id,r)}.groupTuple(by:0)
 filtfa=BLAST_FASTA.out.fasta.filter{m,r,f->!m.control}.map{m,r,f->tuple(m.project_id,f)}.groupTuple(by:0)
 reads_p=PREPARE_QUERIES.out.reads.filter{m,f->!m.control}.map{m,f->tuple(m.project_id,f)}.groupTuple(by:0)
 ustat=MERGE_SEQKIT_STATS.out.stats.filter{p,l,f->l=='FASTQ_unfiltered_statistics'}.map{p,l,f->tuple(p,f)}
 fstat=MERGE_SEQKIT_STATS.out.stats.filter{p,l,f->l=='FASTQ_statistics'}.map{p,l,f->tuple(p,f)}
 mq=MULTIQC_CUSTOMER.out.report.filter{m,f->m.kind=='customer'}.map{m,f->tuple(m.project_id,f)}
 mqd=MULTIQC_CUSTOMER.out.data.filter{m,d->m.kind=='customer'}.map{m,d->tuple(m.project_id,d)}

 // Control/NK/PK reports are copied under PROJECTS/<project>/internal/CONTROLS.
 irawrep=EXCEL_REPORT.out.report.filter{p,mode,kind,x,t->mode=='RAW'&&kind=='internal'}.map{p,m,k,x,t->tuple(p,x,t)}
 ifilrep=EXCEL_REPORT.out.report.filter{p,mode,kind,x,t->mode=='FILTERED'&&kind=='internal'}.map{p,m,k,x,t->tuple(p,x,t)}
 irawhtml=KRONA_HTML.out.html.filter{p,m,k,h->m=='RAW'&&k=='internal'}.map{p,m,k,h->tuple(p,h)}
 ifilhtml=KRONA_HTML.out.html.filter{p,m,k,h->m=='FILTERED'&&k=='internal'}.map{p,m,k,h->tuple(p,h)}
 imq=MULTIQC_CUSTOMER.out.report.filter{m,f->m.kind=='internal'}.map{m,f->tuple(m.project_id,f)}
 imqd=MULTIQC_CUSTOMER.out.data.filter{m,d->m.kind=='internal'}.map{m,d->tuple(m.project_id,d)}
 internalb0=BLAST_CLASSIFY.out.results.filter{m,rawt,filt,assign,hits,ra,fa->m.control}.map{m,rawt,filt,assign,hits,ra,fa->tuple(m.project_id,[rawt,filt,assign,hits])}.groupTuple(by:0).map{p,lists->tuple(p,lists.flatten())}
 internalbqc=BLAST_CLASSIFY.out.qc.filter{m,q->m.control}.map{m,q->tuple(m.project_id,q)}.groupTuple(by:0)
 internalb=internalb0.join(internalbqc,by:0).join(preparation_internal,by:0).map{p,files,qcs,dirs->tuple(p,files+qcs+dirs)}
 rd=README_DOCX.out.docx

 package_input = rawrep
    .join(filrep, by: 0)
    .join(rawhtml, by: 0)
    .join(filhtml, by: 0)
    .join(rawk, by: 0)
    .join(filtk, by: 0)
    .join(rawb, by: 0)
    .join(filtb, by: 0)
    .join(rawfa, by: 0)
    .join(filtfa, by: 0)
    .join(reads_p, by: 0)
    .join(ustat, by: 0)
    .join(fstat, by: 0)
    .join(mq, by: 0)
    .join(mqd, by: 0)
    .join(irawrep, by: 0)
    .join(ifilrep, by: 0)
    .join(irawhtml, by: 0)
    .join(ifilhtml, by: 0)
    .join(imq, by: 0)
    .join(imqd, by: 0)
    .join(internalb, by: 0)
    .join(rd, by: 0)
 analysis_map_for_packaging = PREPARE_ANALYSIS_MAP.out.map.first()
 PACKAGE_PROJECT(package_input,analysis_map_for_packaging)
}
