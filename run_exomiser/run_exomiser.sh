set -eou pipefail
TOOLSPATH=/path/to/exomiser/download/directory/
WORKINGPATH=/path/to/exomiser/yml/file/for/individual
PROPERTIESPATH=/path/to/folder/with/application.properties/file/
STDOUT=$WORKINGPATH/ID.RUNTYPE.stdout
STDERR=$WORKINGPATH/ID.RUNTYPE.stderr ## edit with sample id and exomiser run type name of your choice

module load openjdk/17.0.1
java -jar $TOOLSPATH/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar \
  --analysis $WORKINGPATH/ID_runtype.yml \ ### CHANGE TO NAME OF YML FILE
  --spring.config.location=$PROPERTIESPATH/application.properties \
  > $STDOUT \
  2> $STDERR
