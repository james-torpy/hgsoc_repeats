# load packages:
import(os)
import(glob)

# define starting variables:
project = 'hgsoc_repeats'
expName = 'exp9'

# define directories:
homeDir = '/share/ScratchGeneral/jamtor/'
projectDir = homeDir + '/projects/' + project
resultsDir = projectDir + '/RNA-seq/results'
plotDir = resultsDir + '/R/' + expName + '/plots/DEplots/' + 
descrip + '/'

# fetch filenames:
f_names = glob.glob(plotDir + '/custom3*.pdf')

# fetch samplenames in order:
for l, i in enumerate(file.read(plotDir, + '/s_names.txt')):
	new_name = re.sub('custom3_.*_', 'custom3_' + l + '_', f_names[i])
	os.rename(f_names[i], new_name)