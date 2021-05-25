###################################################################################################
## Parse Input
###################################################################################################
OUT_FOLDER = "results"
SAMPLE_KEY = "example/sampleKey.txt"
TMP_DIR = "~/ad-omics/ricardo/tmp/"

REFERENCE_FASTA = "~/ad-omics/ricardo/Data/1000G_phase1/human_g1k_v37.fasta"
DICT = "~/ad-omics/ricardo/Data/1000G_phase1/human_g1k_v37.dict"
NMASK = "~/ad-omics/ricardo/Data/1000G_phase1/human_g1k_v37.mask.36.fasta.bed"
REF_BUILD = "37"
TOOLS = ["manta","smoove","delly","delly_cnv"]

try:
	f_input_list = pandas.read_csv(SAMPLE_KEY, sep='\t')
	print("Entries in input file:", len(f_input_list.index))
	print("Number of unique samples:", f_input_list['participant_id'].nunique())
except:
	print("Check your input file.")
	raise

participant_id = f_input_list['participant_id'].tolist()
bam_path = f_input_list['bam']

###################################################################################################
## Create folder structure
###################################################################################################
print("Creating folder structure into:", OUT_FOLDER)
os.system("mkdir -p " + OUT_FOLDER)

os.system("mkdir -p " + OUT_FOLDER + "/input")
for label,row in f_input_list.iterrows():
	if os.path.isfile(OUT_FOLDER + "/input/" + row['participant_id'] + ".bam"):
		print("Link to BAM already created: " + row['participant_id'])
	else:
		os.system("ln -f -s " + row['bam'] + " " + OUT_FOLDER + "/input/" + row['participant_id'] + ".bam")


shell.prefix("CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh;")
