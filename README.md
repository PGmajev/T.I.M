**Transcript Intersect & Merge**, or T.I.M for short, is a tool to discover IGRs (InterGenic Regions) in Bacteria. It uses a transcriptome assembly mapped to a reference genome and the annotation of said genome to find features in the transcriptome assembly not found in and not intersecting with annotations in the reference genome.<br />
<br />
**Needed Arguments**:<br />
  --assembly A        The mapped assembly in gff format (filepath from pwd or absolute)<br />
  --genome G          The annotated reference genome in gff3 format (filepath from pwd or absolute)<br />
  --output_path O     Path of the output files (filepath from pwd or absolute)<br />
  --project_prefix N  Project name, all output files will be prefixed with this<br />
<br />
**Optional Arguments**:<br />
  --help, -h          show this help message and exit<br />
  --ss N              Should features only be excluded, if the intersection occurs on their strand? 1= Yes 0=No (Default:0)<br />
  --distance D        Distance of kept feature to annotated features<br />
  --merge M           Should non-intersecting features be merged to the first feature of the merge? 0= Yes 1=No (Default:0)<br />
<br />
Example Command: python3 TIM.py --assembly test_data/Assembly.gff3 --genome test_data/Genome.gff3 --distance 10<br />
<br />
The Input Files have to be gff3 files with the following column order: 'seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes'. While the names are hardcoded, the order of columns may be changed in line 45 of the TIM.py script. <br />
<br />
