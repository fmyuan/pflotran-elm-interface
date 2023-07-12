python3 generateDFN.py
python3 mapdfn2pflotran.py ./output
pflotran -input_prefix cpm_flow
pflotran -input_prefix cpm_transport
