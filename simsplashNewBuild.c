#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

int MAX = 100000;


typedef struct configurationFile{
	char mode[100];
	char name[100];
	char resultsDir[100];
    char starModels[100];
	char starModelsDir[100];
    char modelType[100];
	char rungeKuttaType[100];
	char range[100];
	char makePlot[100];
	char makeHR[100];
    double initialSMA;
    double planetMass;
	double starMass;
	double starMetallicity;
	

} configuration;


typedef struct starEvolution{
	double massStar;
	double envelopemassStar;
	double radiusStar; 
	double luminosityStar; 
	double masslossStar; 
	double phaseStar; 
	double ageStar; 
	double rotationStar; 
	double coremassStar; 
	double helium4Star; 
	double hydrogen1Star;
	double newSemiMajorAxis;
	struct starEvolution *next;
} starEvolution;

configuration *readConfig();
starEvolution *createEvolutionData(starEvolution *evolveSemiMajorAxis);
int semiMajorAxisEvolution(configuration *simulation);
int cleanStellarModel(configuration *simulation);
int evolveSystems(configuration *simulation);

int main(){	
	
	if( access( "./configuration", F_OK ) == -1 ) {
		printf("ERROR: The configuration file is missing. Please place the configuration in working directory\n");
		return 0;	}
	mkdir("./Results", 0777);
	struct configurationFile *simulation = NULL;	
	simulation = readConfig();
	
	if(strcmp((simulation->mode), "semiMajorAxisEvolution") == 0){
		semiMajorAxisEvolution(simulation);}
	if(strcmp((simulation->mode), "singleEpochPopulation") == 0){
		printf("WARNING: The single epoch mode is currently under refurbishment. Please try using older version of simsplash or wait for new version\n");}
	if(strcmp((simulation->mode), "multiEpochPopulation") == 0){
		printf("WARNING: The multi epoch mode is currently under refurbishment. Please try using older version of simsplash or wait for new version\n");}
		
	return 0;
}

configuration *readConfig(){	
	configuration *simulation = malloc(sizeof(struct configurationFile));
	char *configText = (char *)malloc(sizeof(char));
	char *resultsDir = (char *)malloc(sizeof(char));
	char *starModelsDir = (char *)malloc(sizeof(char));
	char *field = (char *)malloc(sizeof(char));
	char *value = (char *)malloc(sizeof(char));
	FILE *INPUT = fopen("./configuration", "r");
	//To do: check whether I can just use grep for this 
	while(!feof(INPUT)){		
		fscanf(INPUT,"%[^\n]%*c\n", &*configText);				
		if(configText[0] != '#'){
			sscanf(configText, "%s = %s", &*field, &*value);
			if (strcmp(field,"Mode") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->mode);}
			if (strcmp(field,"Name") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->name);}
			if (strcmp(field,"starModels") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->starModels);}				
			if (strcmp(field,"modelType") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->modelType);}				
			if (strcmp(field,"rungeKuttaType") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->rungeKuttaType);}	
			if (strcmp(field,"rangeSMA") == 0){	
				sscanf(configText, "%s = %s", &*field, simulation->range);}				
			if (strcmp(field,"initialSMA") == 0){	
				sscanf(configText, "%s = %lf", &*field, &simulation->initialSMA);}			
			if (strcmp(field,"starMass") == 0){	
				sscanf(configText, "%s = %lf", &*field, &simulation->starMass);}
			if (strcmp(field,"planetMass") == 0){	
				sscanf(configText, "%s = %lf", &*field, &simulation->planetMass);}
			if (strcmp(field,"starMetallicity") == 0){	
				sscanf(configText, "%s = %lf", &*field, &simulation->starMetallicity);}
			if (strcmp(field,"makePlot") == 0){	
				sscanf(configText, "%s = %s", &*field, &simulation->makePlot);}	
			if (strcmp(field,"makeHR") == 0){	
				sscanf(configText, "%s = %s", &*field, &simulation->makeHR);}				
		}		
	}
	sprintf(resultsDir, "./Results/%s", simulation->name);
	sprintf(starModelsDir, "../Stellar_Database/%s/%gz_stars/%.1f_%g", simulation->starModels, simulation->starMetallicity, simulation->starMass, simulation->starMetallicity);

	sscanf(starModelsDir, "%s", simulation->starModelsDir);
	sscanf(resultsDir, "%s", simulation->resultsDir);
	return simulation;
	
}

int semiMajorAxisEvolution(configuration *simulation){
	char error[] = "\n\nERROR: There is a problem with the simulation initialisation. \nPlease check the configuration file";
	int status = 1;
	printf("******You have chosen the Semi Major Axis Evolution Mode of SIMSPLASH******\n\nChecking Stellar Models....\n\n");	
	
	status = cleanStellarModel(simulation);
	if (status == 0){
		return 0;}
		
	if(strcmp(simulation->modelType, "single") == 0){	
	evolveSystems(simulation);
	}
	else if(strcmp(simulation->modelType, "range") == 0){	
	//Add range of semi major axis here, for loop calling evolvesystems, where the struct->initiala changes on every iteration
	evolveSystems(simulation);
	}
	else{
		printf("%s\n", error);
		return 0;
	}
	
	
}

int cleanStellarModel(configuration *simulation){
	char stellarModelError[] = "*\n*\n*\n*\n*\nATTENTION:: ERROR: The specified stellar models do not exist.\nCheck configuration file\n";
	char *cleanHistoryFile = (char*)malloc(sizeof(char)*MAX);
	char *extractHistoryFields = (char*)malloc(sizeof(char)*MAX);
	
	if (strcmp(simulation->starModels,"MIST_Rotating")==0 || strcmp(simulation->starModels,"MIST_NonRotating")==0){
		if ((simulation->starMetallicity == 0.005) || (simulation->starMetallicity == 0.008) || (simulation->starMetallicity == 0.014) || (simulation->starMetallicity == 0.025) || (simulation->starMetallicity == 0.04)){
			mkdir(simulation->resultsDir, 0777);
			printf("New results directory created as: %s\n", simulation->resultsDir);
		}
		else{ printf("%s", stellarModelError);
			return 0;
		}	
	}
	else if (strcmp(simulation->starModels,"MESA_Rotating")==0 || strcmp(simulation->starModels,"MESA_NonRotating")==0){
		if ((simulation->starMetallicity == 0.01) || (simulation->starMetallicity == 0.02)){
			mkdir(simulation->resultsDir, 0777);
			printf("New results directory created as: %s\n", simulation->resultsDir);
		}
		else{ printf("%s", stellarModelError);
			return 0;
		}
	}	
	else{ printf("%s", stellarModelError);
			return 0;
		}

	sprintf(cleanHistoryFile, "sed 's/#//' %s/history.data > %s/starData.data; sed -i -n -e '/star_age/,$p' %s/starData.data", simulation->starModelsDir, simulation->starModelsDir, simulation->starModelsDir);
	system(cleanHistoryFile);
	if (strcmp(simulation->starModels,"MIST_Rotating")==0 || strcmp(simulation->starModels,"MIST_NonRotating")==0){
		sprintf(extractHistoryFields, "awk ' {if(NR==1) for(i=1;i<=NF;i++) { ix[$i] = i} if(NR > 1) print $ix[c1], $ix[c2], $ix[c3], $ix[c4], $ix[c5], $ix[c6], $ix[c7], $ix[c8]}' c1=star_mass c2=star_age c3=star_mdot c4=log_L c5=log_R c6=he_core_mass c7=surf_avg_omega c8=phase %s/starData.data > %s/starData.txt", simulation->starModelsDir, simulation->starModelsDir );	
	}
	if (strcmp(simulation->starModels,"MESA_Rotating")==0 || strcmp(simulation->starModels,"MESA_NonRotating")==0){
		sprintf(extractHistoryFields, "awk ' {if(NR==1) for(i=1;i<=NF;i++) { ix[$i] = i} if(NR > 1) print $ix[c1], $ix[c2], $ix[c3], $ix[c4], $ix[c5], $ix[c6], $ix[c7], $ix[c8], $ix[c9]}' c1=star_mass c2=star_age c3=star_mdot c4=log_L c5=log_R c6=he_core_mass c7=surf_avg_omega c8=center_h1 c9=centerh4 %s/starData.data > %s/starData.txt", simulation->starModelsDir, simulation->starModelsDir );
	}
	
	system(extractHistoryFields);
	return 1;
}

int evolveSystems(configuration *simulation){
	char *starData = (char *)malloc(sizeof(char)*MAX);
	starEvolution *evolveSemiMajorAxis = malloc(sizeof(struct starEvolution));
	sprintf(starData, "%s/starData.txt", simulation->starModelsDir);
	FILE *INPUT = fopen(starData, "r");
	int phaseCounter = 0;
	
	evolveSemiMajorAxis->newSemiMajorAxis = simulation->initialSMA;
	
	while(!feof(INPUT)){
		fscanf(INPUT, "%lf %lf %lf %lf %lf %lf %lf %lf", &evolveSemiMajorAxis->massStar, &evolveSemiMajorAxis->ageStar, &evolveSemiMajorAxis->masslossStar, &evolveSemiMajorAxis->luminosityStar, &evolveSemiMajorAxis->radiusStar, &evolveSemiMajorAxis->coremassStar, &evolveSemiMajorAxis->rotationStar, &evolveSemiMajorAxis->phaseStar);
		evolveSemiMajorAxis->phaseStar = round(evolveSemiMajorAxis->phaseStar);				
		evolveSemiMajorAxis->radiusStar = pow(10.0, evolveSemiMajorAxis->radiusStar);
		evolveSemiMajorAxis->luminosityStar = pow(10.0, evolveSemiMajorAxis->luminosityStar);			
		evolveSemiMajorAxis->envelopemassStar = evolveSemiMajorAxis->massStar-evolveSemiMajorAxis->coremassStar;
		

		
	}
		
	return 1;
}

starEvolution *createEvolutionData(starEvolution *evolveSemiMajorAxis){
	
}
















