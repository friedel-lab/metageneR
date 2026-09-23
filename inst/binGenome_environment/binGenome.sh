#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='java:rm' # Rscript
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

INSTALL_DIR="$SCRIPT_FOLDER/R/"
SCRIPT_PATH="binGenome.R"

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'bedgraph' '' "one ore more bedgraph files separated by ','" 'i'
DEFINE_string 'bedgraphPos' '' "one ore more bedgraph files for positive strand separated by ','" ''
DEFINE_string 'bedgraphNeg' '' "one ore more bedgraph files for negative strand separated by ','" ''
DEFINE_string 'annotation' '' "one ore more region annotation files separated by ','; (see writeGRangesToBed() in R/binGenome.lib.R  for format" 'a'
DEFINE_string 'bins' '' "[optional] one ore more number of bins separated by ','" 'b'
DEFINE_string 'quantiles' '' "[optional] determines the position at which expression exceeds specific quantiles; separated by ','" 'q'
DEFINE_string 'outputDir' '' 'path to output folder; files will be named automatically based on the used parameters' 'o'
DEFINE_string 'bedgraphNames' '' "[optional] sample names for output filenames separated by ','" ''
DEFINE_string 'annotationNames' '' "[optional] annotation names for output filenames separated by ','" ''
DEFINE_integer 'cores' '1' '[optional] number of cores to use in parallel' 'c'
DEFINE_boolean 'normalize' 'true' '[optional] write in addition a per-gene normalized version of the data' 'n'
DEFINE_string 'fixedBinSizeUpstream' '' "[optional] can be used to create fixed bins upstream; format: 'binsize:binnumber'" ''
DEFINE_string 'fixedBinSizeDownstream' '' "[optional] can be used to create fixed bins downstream; format: 'binsize:binnumber'" ''
DEFINE_string 'tmpDir' '' "[optional] path to tmp folder" ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_bedgraph" ] && [ -z "$FLAGS_bedgraphPos" ] && [ -z "$FLAGS_bedgraphNeg" ]; then
	echoError "Parameter --bedgraph or --bedgraphPos and --bedgraphNeg must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_annotation" ]; then
	echoError "Parameter --annotation must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_bins" ] && [ -z "$FLAGS_quantiles" ] ; then
	echoError "Parameter --bins or --quantiles must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_outputDir" ]; then
	echoError "Parameter --outputDir must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

### all further checks are done by the Rscript ###

# all parameters seem to be ok, run the script!
RUN2END=$(getTmpFile "binGenome")
MEMORY=$(getMemoryForJava $FLAGS_threads 5012)
COMMAND="java -XX:-UsePerfData -Djava.io.tmpdir=/mnt/biocluster/studtemp/tabanli/genome_binner $MEMORY -jar '$SCRIPT_FOLDER/genomeBinner.jar'"
#/mnt/biocluster/studtemp/tabanli/genome_binner
#/home/t/tabanli/tmp/genomeBinner_tmp
#/tmp

# build the command
for PARAM in $__flags_longNames; do
	# ignore that parameter since it is only for the module
	if [ "$PARAM" == "debug" ]; then
		continue
	fi

	V_NAME='FLAGS_'"$PARAM"
	VALUE=${!V_NAME}
	V_NAME='__flags_'"$PARAM"'_type'
	TYPE=${!V_NAME}
	
	# test if value for that parameter was set
	if [ ! -z "$VALUE" ]; then
		# test if flag
		if [ $TYPE -eq 1 ]; then
			# test if flag is enabled
			if  [ $VALUE -eq 0 ]; then
				COMMAND="$COMMAND --$PARAM"
			fi
		else 
			COMMAND="$COMMAND --$PARAM '$VALUE'"
		fi
	fi
done
echo $COMMAND
# execute the command
#LOG="$FLAGS_outputDir/Rscript.log" # OLD R call
mkdir -p "$FLAGS_outputDir"
LOG="$FLAGS_outputDir/GenomeBinner_"$(date +%s_%N)".log"
echo "LOG path: $LOG"
if [ -e "$LOG" ]; then
	rm "$LOG"
fi
executeCommand "$COMMAND" "$LOG" "java genomeBinner" 
# test, if the script run to its end

# if, we come till here, all should be ok.
exit $EXIT_OK
