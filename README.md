**Needed Arguments**:<br />
  --assembly A        The mapped assembly in gff format (filepath from pwd or absolute)<br />
  --genome G          The annotated reference genome in gff3 format (filepath from pwd or absolute)<br />
  --output_path O     Path of the output files (filepath from pwd or absolute)<br />
  --project_prefix N  Project name, all output files will be prefixed with this<br />
<br />
**Optional Arguments**:<br />
  -h, --help          show this help message and exit<br />
  --ss N              Should features only be excluded, if the intersection occurs on their strand? 1= Yes 0=No (Default:0)<br />
  --distance D        Distance of kept feature to annotated features<br />
  --merge M           Should non-intersecting features be merged to the first feature of the merge? 0= Yes 1=No (Default:0)<br />

