
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "paquetReading.h"
#include "interpreter.h"
#include "image.h"
#include "data.h"

// variables globales
char *g_dataModelPath;
char *g_dataFilePath;

const corr CHANNEL_RESOLUTION[17];

/*  Check and exit
 */
void Check(int status, int ligne, char* fichier){
	if (status == -3) {
		printf("Erreur East %d : fichier modèle de données .eas non reconnu \n",status);
		exit(1);
	}
	if (status != 1) {
		printf("Erreur East %d : ligne %d du fichier %s \n",status,ligne,fichier);
		exit(1);
	}
}

/* - Return the value of the width of the widest image (useful for creating the full disk)
 *
 *  - input : array of images
 *  			 number of images in the array
 */
int MaxTab(Image** images,int nbrImages){
	int i, MAX;
	MAX = 0;
	for (i=0; i<nbrImages; i++)
		if(images[i]->w > images[MAX]->w ) MAX = i;
	return images[MAX]->w;
}

/* - Return the number of file in a directory
 *
 *  - input : name of directory (within /test_best/bin/workspace/Project1/.. )
 */
int GetNumberOfFile(char directory[50]){
	int k = 0, NB_FILE;
	char tmp[1][4];
	char chaineCommand[300] = "cd ";
	strcat(chaineCommand,directory);
	strcat(chaineCommand," ; ls -1 | wc -l");

	//get the number of files in the directory FCI_files
	FILE *command = popen(chaineCommand, "r");
	while( fgets(tmp[k], 4, command) != NULL ){
		k++;
	}
	pclose(command);
	NB_FILE = strtol(tmp[0],NULL,10);

	return NB_FILE;
}

/* - Stock as string the full path of all files of a directory
 *
 *  - input : an array of string = where all path will be stored
 *  			 an int = number of files in the directory
 *  			 a string = absolute path of where the files are
 */
void GetNameFiles(char** path, int NB_FILE, char *fullPath){

	int k = 0, i, taille;
	char chaineCommand1[300] = "cd ";
	strcat(chaineCommand1,fullPath);
	char chaineCommand2[300] = "; readlink -f $(ls ";
	strcat(chaineCommand2,fullPath);
	strcat(chaineCommand2,")");
	strcat(chaineCommand1,chaineCommand2);

	//get the path of our binary files in an array of string
	FILE *command = popen(chaineCommand1, "r");
	k = 0;
	while( fgets(path[k], 500, command) != NULL ){
		k++;
	}
	pclose(command);
	for(i=0 ; i < NB_FILE ; i++) { // erase \n at the end of our string
		taille = strlen(path[i]);
		path[i][(taille - 1)] = '\0';
	}

}

/* - Get the Scan angles and write them in a file
 *
 *  - input : an int = the swath
 */
void WriteScanAngles(int a_swath,int a_RC){

	SCAaxis *NSaxis, *EWaxis;
	int **fineTime, k, nbrSCA = 0, i, j, status;
	char nameFichierScan[50] = "ScanAngles/RCxx_sw";
	char tmpRC[3],tmpsw[3];
	sprintf(tmpRC, "%d", a_RC);
	if (LengthInteger(a_RC) == 1){
		tmpRC[1] = tmpRC[0];
		tmpRC[0] = '0';
		tmpRC[2] = '\0';
	}
	nameFichierScan[13] = tmpRC[0];nameFichierScan[14] = tmpRC[1];
	sprintf(tmpsw, "%d", a_swath);
	if (LengthInteger(a_swath) == 1){
		tmpsw[1] = tmpsw[0];
		tmpsw[0] = '0';
		tmpsw[2] = '\0';
	}
	strcat(nameFichierScan,tmpsw);
	strcat(nameFichierScan,".dat");

	NSaxis = malloc(10000 * sizeof(*NSaxis));
	EWaxis = malloc(10000 * sizeof(*EWaxis));
	fineTime = malloc(10000 * sizeof(*fineTime));
	for (k=0 ; k < 10000 ; k++)
		fineTime[k] = malloc(50 * sizeof(**fineTime));

	status = GetScanAngles(NSaxis,EWaxis,fineTime,&nbrSCA);
	Check(status,__LINE__-1,__FILE__);
	FILE* fichierEcrire = NULL;
	fichierEcrire = fopen(nameFichierScan,"w");
	if (fichierEcrire != NULL) {
		fprintf(fichierEcrire,"    NS       EW\n");
		for (i=0;i<nbrSCA;i++) { // lines
			for (j=0;j<50;j++) { // 50 angles by packets
				fprintf(fichierEcrire,"%f ",NSaxis[i].angle[j]);
				fprintf(fichierEcrire,"%f \n",EWaxis[i].angle[j]);
			}
		}
		fclose(fichierEcrire);
		printf("fichier \"%s\" écrit !\n",nameFichierScan);
	}

	for(i=0 ; i < 10000 ; i++){
		free(fineTime[i]);
	}
	free(fineTime);
	free(NSaxis);
	free(EWaxis);

}

/* Build image of the LI
 *
 * input : an int = file number of the image to be extracted = order in which it appears in the directory
 * 		   an int = number of the Optical Channel , 0 if we want full image
 */
void GetImageLI(int a_RC, int a_OC){

	//variables
	int status;
	char tmp[10];
	char nameImageOC[100] = "imagesLI/RCxx_OCx.bmp";

	//create name of images according to directory = file asked in argument
	sprintf(tmp,"%d",a_RC);
	if (LengthInteger(a_RC) == 1){
		nameImageOC[11] = '0';
		nameImageOC[12] = tmp[0];
	} else {
		nameImageOC[11] = tmp[0];
		nameImageOC[12] = tmp[1];
	}

	//create name of images according to Optical Channel number
	sprintf(tmp,"%d",a_OC);
	nameImageOC[16] = tmp[0];

	if ( a_OC == 0 ) {
		status = GetFullOC(nameImageOC);
		Check(status,__LINE__-1,__FILE__);
	} else {
		status = GetOC(a_OC, nameImageOC);
		Check(status,__LINE__-1,__FILE__);
	}


}

/* Build image bitmap of a swath of the FCI by calling GetChannel() which get all value of pixel in an array and calling CreateImage()
 * the other operation are manage name of the file and his size
 *
 * input : an int = swath's number
 * 		   char* = channel's name
 * 		   int = resolution
 */
void GetImageFCI(int a_RC, int a_swath, char *a_nameChannel, int a_resolution){

	int status, i, k, size, nbrPixelByColumn = 0, nbrPixelByRow = 0;
	char tmpRC[4],tmpsw[4], tmpres[4];
	sprintf(tmpRC, "%d", a_RC);
	if (LengthInteger(a_RC) == 1){
		tmpRC[0] = '0';
		tmpRC[1] = tmpRC[0];
	}
	sprintf(tmpsw, "%d", a_swath);
	if (LengthInteger(a_swath) == 1){
		tmpsw[1] = tmpsw[0];
		tmpsw[0] = '0';
		tmpsw[2] = '\0';
	}
	sprintf(tmpres, "%d",  a_resolution);


	if (strcmp(a_nameChannel,"all") == 0){

		int ***channel;
		nbrPixelByColumn = 448;

		channel = malloc(17 * sizeof(*channel));
		if (channel == NULL)
			printf("erreur d'allocation : ligne %d du fichier %s \n",__LINE__,__FILE__);
		for(k=0 ; k < 17 ; k++){
			channel[k] = malloc(nbrPixelByColumn * sizeof(**channel));
			if (channel[k] == NULL)
				printf("erreur d'allocation : ligne %d du fichier %s \n",__LINE__,__FILE__);
		}
		for(k=0 ; k < 17 ; k++){
			for(i=0 ; i < nbrPixelByColumn ; i++){
				channel[k][i] = malloc(100000 * sizeof(***channel));
				if (channel[k][i] == NULL)
					printf("erreur d'allocation : ligne %d du fichier %s \n",__LINE__,__FILE__);
			}
		}
		int taille1, taille2, taille3;
		status = GetAllChannel(channel,a_resolution, &taille1, &taille2, &taille3);
		Check(status,__LINE__-1,__FILE__);

		for (k=0;k<17;k++){

			char nameImage[100] = "imagesFCI/RCxx_";
			nameImage[12] = tmpRC[0]; nameImage[13] = tmpRC[1];
			strcat(nameImage,CHANNEL_RESOLUTION[k].channel);
			strcat(nameImage,"_sw");
			strcat(nameImage,tmpsw);
			strcat(nameImage,"_res");
			strcat(nameImage,tmpres);
			strcat(nameImage,".bmp");

			size = CHANNEL_RESOLUTION[k].size;
			nbrPixelByColumn = size;

			if (size == 112)
				nbrPixelByRow = taille3;
			else if (size == 224)
				nbrPixelByRow = taille2;
			else if (size == 448)
				nbrPixelByRow = taille1;

			if (CHANNEL_RESOLUTION[k].channel[0]=='V') {
				if (a_swath%2 == 0)
					CreateImage(channel[k], nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "noflip");
				else
					CreateImage(channel[k], nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip horizontal");
			}else{
				if (a_swath%2 == 0)
					CreateImage(channel[k], nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip vertical");
				else
					CreateImage(channel[k], nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip horizontal vertical");
			}
		}

		for(k=0 ; k < 17 ; k++){
			for(i=0 ; i < nbrPixelByColumn ; i++){
				free(channel[k][i]);
			}
			free(channel[k]);
		}
		free(channel);

	} else {

		int **channel;

		for (k=0;k<17;k++){

			if (strcmp(CHANNEL_RESOLUTION[k].channel,a_nameChannel) == 0)
				size = CHANNEL_RESOLUTION[k].size;
		}
		// create 2D array of int (where image is stocked) according to his size = argument resolution
		nbrPixelByColumn = size;

		channel = malloc(nbrPixelByColumn * sizeof(*channel));
		for(i=0 ; i < nbrPixelByColumn ; i++)
			channel[i] = malloc(100000 * sizeof(**channel));

		status = GetChannel(channel, a_nameChannel, &nbrPixelByRow, a_resolution); //get value of pixel of the channel
		Check(status,__LINE__-1,__FILE__);


		// create swath image and save it
		char nameImage[100] = "imagesFCI/RCxx_";
		nameImage[12] = tmpRC[0]; nameImage[13] = tmpRC[1];
		strcat(nameImage,a_nameChannel);
		strcat(nameImage,"_sw");
		strcat(nameImage,tmpsw);
		strcat(nameImage,"_res");
		strcat(nameImage,tmpres);
		strcat(nameImage,".bmp");

		if (a_nameChannel[0]=='V') {
			if (a_swath%2 == 0)
				CreateImage(channel, nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "noflip");
			else
				CreateImage(channel, nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip horizontal");
		}else{
			if (a_swath%2 == 0)
				CreateImage(channel, nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip vertical");
			else
				CreateImage(channel, nameImage, nbrPixelByRow/a_resolution, nbrPixelByColumn/a_resolution, "flip horizontal vertical");
		}

		for(i=0 ; i < nbrPixelByColumn ; i++){
			free(channel[i]);
		}
		free(channel);

	}

}

/* Create and Assemble each swath images one above another
 *
 * input :  char* = name of the channel
 * 			int = number of the 1st swath
 * 			int = number of the last swath
 * 			int = resolution of swath
 */
void AssembleImages(int a_RC,int a_posImage1, int a_posImage2, char* a_nameChannel,int a_resolution){

	int i, swath , nbrPixelByColumn = 0, nbrPixelByRow, height, width, Uturn;
	int nbrImages = a_posImage2 - a_posImage1;
	char *flipped;
	Image **images;
	images = malloc(nbrImages* sizeof(*images));
	Image* imageAssemblee = NULL;
	char tmpres[4],tmpRC[4];
	sprintf(tmpres, "%d",  a_resolution);
	sprintf(tmpRC, "%d", a_RC);
	if (LengthInteger(a_RC) == 1){
		tmpRC[1] = tmpRC[0];
		tmpRC[0] = '0';
	}

	printf("%d images\n",nbrImages);

	int NB_FILE = GetNumberOfFile(g_dataFilePath); //stock name of datafile in an array
	char**	path = malloc(NB_FILE * sizeof(*path));
	for(i=0 ; i < NB_FILE ; i++)
		path[i] = malloc(500 * sizeof(**path));
	GetNameFiles(&*path,NB_FILE,g_dataFilePath);

	for (swath = a_posImage1;swath<a_posImage2;swath++){
		char tmpsw[5];
		sprintf(tmpsw, "%d", swath);
		if (LengthInteger(swath) == 1){
			tmpsw[1] = tmpsw[0];
			tmpsw[0] = '0';
			tmpsw[2] = '\0';
		}

		char nameImage[100] = "imagesFCI/RCxx_";
		nameImage[12] = tmpRC[0]; nameImage[13] = tmpRC[1];
		strcat(nameImage,a_nameChannel);
		strcat(nameImage,"_sw");
		strcat(nameImage,tmpsw);
		strcat(nameImage,"_res");
		strcat(nameImage,tmpres);
		strcat(nameImage,".bmp");
		g_dataFile = path[swath];

		images[swath-a_posImage1] = Charger(nameImage);
		if (images[swath-a_posImage1] == NULL) {
			GetImageFCI(a_RC, swath, a_nameChannel,a_resolution);
		}

	}
	int k,size;
	for (k=0;k<17;k++){
		if (strcmp(CHANNEL_RESOLUTION[k].channel,a_nameChannel) == 0)
			size = CHANNEL_RESOLUTION[k].size;
	}
	if (size == 112) // taille du Uturn à enlever dans le swath
		Uturn = 2089/(4*a_resolution); // empiriquement ?
	else if (size == 224)
		Uturn = 2089/(2*a_resolution);
	else if (size == 448)
		Uturn = 2089/(1*a_resolution);

	//initialize new image
	nbrPixelByRow = MaxTab(images,nbrImages)+0.1*MaxTab(images,nbrImages);//-Uturn;
	for (i=0;i<nbrImages;i++)
		nbrPixelByColumn = nbrPixelByColumn+images[i]->h;

	imageAssemblee = NouvelleImage(nbrPixelByRow,nbrPixelByColumn);

	width = imageAssemblee->w/2-((images[nbrImages-1]->w-Uturn)/2); // where to start
	height = 0;

	if (a_posImage2 % 2 == 0)
		flipped = "flip";
	else
		flipped = "noflip";

	AddImageNoUturn(images[nbrImages-1],imageAssemblee,width,height,flipped,Uturn);

	for (i = nbrImages-2;i>-1;i--){

		width = imageAssemblee->w/2-((images[i]->w-Uturn)/2);
		height = height+(images[i+1]->h);

		if (a_posImage1 % 2 == 0) {
			if (i % 2 == 0)
				flipped = "noflip";//		i = Uturn : image.w
			else
				flipped = "flip";//			i = 0 : image.w-Uturn
		} else {
			if (i % 2 == 0)
				flipped = "flip";
			else
				flipped = "noflip";
		}

		AddImageNoUturn(images[i],imageAssemblee,width,height,flipped,Uturn); //add an image at position x = width without the Uturn
	}

	char nameFullImage[100] = "imagesFCI/RCxx_FullDisk_";
	nameFullImage[12] = tmpRC[0]; nameFullImage[13] = tmpRC[1];
	strcat(nameFullImage,a_nameChannel);
	strcat(nameFullImage,"_res");
	strcat(nameFullImage,tmpres);
	strcat(nameFullImage,".bmp");

	SaveImage(imageAssemblee,nameFullImage);
	printf("\"%s\" sauvegardée \n",nameFullImage);
	DeleteImage(imageAssemblee);

	for(i=0 ; i < nbrImages ; i++){
		free(images[i]);
	}
	free(images);
}

/* Assemble swath images of all Channel one above another
 *
 * input :  char* = name of the channel
 * 			int = number of the 1st swath
 * 			int = number of the last swath
 * 			int = resolution of swath
 */
void AssembleAllImages(int a_RC,int a_posImage1, int a_posImage2, int a_resolution){

	int k, i, swath , nbrPixelByColumn = 0, nbrPixelByRow = 0, height, width, Uturn,size ;
	int nbrImages = a_posImage2 - a_posImage1;
	char *flipped;
	Image* imageAssemblee = NULL;
	Image ***images;
	images = malloc(17 * sizeof(*images));
	for(i=0 ; i < 17 ; i++)
		images[i] = malloc(nbrImages * sizeof(**images));

	char tmpres[4],tmpRC[4];
	sprintf(tmpres, "%d",  a_resolution);
	sprintf(tmpRC, "%d", a_RC);
	if (LengthInteger(a_RC) == 1){
		tmpRC[0] = '0';
		tmpRC[1] = tmpRC[0];
	}

	int NB_FILE = GetNumberOfFile(g_dataFilePath); //stock path of datafile in array
	char**	path = malloc(NB_FILE * sizeof(*path));
	for(i=0 ; i < NB_FILE ; i++)
		path[i] = malloc(500 * sizeof(**path));
	GetNameFiles(&*path,NB_FILE,g_dataFilePath);

	for (swath = a_posImage1;swath<a_posImage2;swath++){
		g_dataFile = path[swath];

		GetImageFCI(a_RC,swath, "all" ,a_resolution);

		for (k=0;k<17;k++){
			char tmpsw[4];
			sprintf(tmpsw, "%d", swath);
			if (LengthInteger(swath) == 1){
				tmpsw[1] = tmpsw[0];
				tmpsw[0] = '0';
				tmpsw[2] = '\0';
			}
			char nameImage[100] = "imagesFCI/RCxx_";
			nameImage[12] = tmpRC[0]; nameImage[13] = tmpRC[1];
			strcat(nameImage,CHANNEL_RESOLUTION[k].channel);
			strcat(nameImage,"_sw");
			strcat(nameImage,tmpsw);
			strcat(nameImage,"_res");
			strcat(nameImage,tmpres);
			strcat(nameImage,".bmp");
			if (images[k][swath-a_posImage1]==NULL)
				printf("image non chargée !\n");
			images[k][swath-a_posImage1] = Charger(nameImage);

		}
	}

	for (k=0;k<17;k++){

		size = CHANNEL_RESOLUTION[k].size;
		if (size == 112)
			Uturn = 2089/(4*a_resolution);//empiriquement
		else if (size == 224)
			Uturn = 2089/(2*a_resolution);
		else if (size == 448)
			Uturn = 2089/(1*a_resolution);

		nbrPixelByRow = MaxTab(images[k],nbrImages)-Uturn;//+0.1*MaxTab(images,nbrImages);
		nbrPixelByColumn = 0;
		for (i=0;i<nbrImages;i++)
			nbrPixelByColumn = nbrPixelByColumn+images[k][i]->h;
		imageAssemblee = NouvelleImage(nbrPixelByRow,nbrPixelByColumn);

		width = imageAssemblee->w/2-((images[k][nbrImages-1]->w-Uturn)/2); // where to start
		height = 0;

		if (a_posImage2 % 2 == 0)
			flipped = "flip";
		else
			flipped = "noflip";

		AddImageNoUturn(images[k][nbrImages-1],imageAssemblee,width,height,flipped,Uturn);

		for (i = nbrImages-2;i>-1;i--){
			if (a_posImage1 % 2 == 0) {
				if (i % 2 == 0)
					flipped = "noflip";//		i = Uturn : image.w
				else
					flipped = "flip";//	i = 0 : image.w-Uturn
			} else {
				if (i % 2 == 0)
					flipped = "flip";
				else
					flipped = "noflip";
			}

			width = imageAssemblee->w/2-((images[k][i]->w-Uturn)/2);
			height = height+(images[k][i+1]->h);
			AddImageNoUturn(images[k][i],imageAssemblee,width,height,flipped,Uturn);

		}

		char nameFullImage[100] = "imagesFCI/RCxx_FullDisk_";
		nameFullImage[12] = tmpRC[0]; nameFullImage[13] = tmpRC[1];
		strcat(nameFullImage,CHANNEL_RESOLUTION[k].channel);
		strcat(nameFullImage,"_res");
		strcat(nameFullImage,tmpres);
		strcat(nameFullImage,".bmp");

		SaveImage(imageAssemblee,nameFullImage);
		printf("\"%s\" sauvegardée \n",nameFullImage);
		DeleteImage(imageAssemblee);

	}

	for(i=0 ; i < nbrImages ; i++)
		free(images[i]);

	free(images);
}

/* Teste si les parametres rentrés par l'utilisateur sont corrects
 *
 */
int TestParam(char* a_nameChannel, int a_resolution) {
	int testNameChannel = 0, k;

	if (strcmp(a_nameChannel,"all") == 0)
		testNameChannel = 1;

	for (k=0;k<17;k++){
		if (strcmp(CHANNEL_RESOLUTION[k].channel,a_nameChannel) == 0)
			testNameChannel = 1;
	}
	if (testNameChannel == 0 ) {
		printf("Erreur : nom du canal incorrect !\n");
		printf("Choisir parmi : VIS_1, VIS_2, VIS_3, VIS_4, VIS_5, NIR_1, NIR_2, NIR_3, IR1_1, IR1_2, IR2_1, IR2_2, IR3_1, IR3_2, IR3_3");
		return -1;
	}
	if (a_resolution <= 0 ) {
		printf("Erreur : resolution <= 0 impossible !\n");
		return -1;
	}

	return 1;
}

/* Teste si le fichier rentré par l'utilisateur est reconnu ou non
 *
 */
int TestFile(char* a_file){

	if (fopen(a_file, "r+") == NULL){
		printf("Erreur : fichier de données inconnu\n");
		return -1;
	}
	return 1;
}

int main(int argc, char *argv[]) {

	char tmp2[300];
	char tmp1[300];
	int swath, resolution, RC ;
	char* nameChannel = malloc(sizeof(char) * 10);
	int choixMenu = 0, status, quitter = 1,numeroBloc;

	g_dataFile = malloc(sizeof(char) * 500);
	g_dataFilePath = malloc(sizeof(char) * 500);

	while (quitter) {

		time_t depart, arrivee;
		time(&depart);

		g_dataModelPath = argv[1];

		status = start_i(FRENCH);		// Initialisation of services EAST
		Check(status,__LINE__-1,__FILE__);

		parameterMOF parameterMOF;
		parameterMOF.UTC0 = 1450742399.9023139;
		parameterMOF.gradient = 1;
		parameterMOF.OBT0 = 23932900.402313948;
		parameterMOF.offset = 0;

		printf("\n======================    Menu    ========================\n\n");
		printf("[1] FCI   [2] LI   [3] Header   [4] Plateforme   [5] Quitter\n\n");
		scanf("%d", &choixMenu);
		switch (choixMenu)  {
		case 1:

			strcpy(tmp2,g_dataModelPath); // create value of the data model and the data files
			strcat(tmp2,"MODEL_FCI.eas");
			g_dataModel = tmp2;

			status = select_DDR_i(g_dataModel);
			Check(status,__LINE__-1,__FILE__);

			int choixFCI = 0;
			printf("=============================    FCI   =============================\n\n");
			printf("[1] Image swath  [2] Full Disk  [3] Angles de scan  [4] RC start time \n\n");
			scanf("%d", &choixFCI);
			if ((choixFCI != 0) && (choixFCI != 1) && (choixFCI != 2) && (choixFCI != 3) && (choixFCI != 4)){
				printf("Vous n'avez pas rentré un nombre correct.\n");
				break;
			}

			switch (choixFCI)  {
			case 1:
				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				if (TestFile(g_dataFile) != 1)
					break;
				printf("\nNuméro du Repeat Cycle : ");
				scanf("%d", &RC);
				printf("\nNuméro du swath : ");
				scanf("%d", &swath);
				printf("\nNom du canal : ");
				scanf("%s", nameChannel);
				printf("\nRésolution de l'image : ");
				scanf("%d", &resolution);
				printf("\n");
				if (TestParam(nameChannel,resolution) == 1){
					printf("... décodage en cours ...\n");
					GetImageFCI(RC,swath,nameChannel, resolution);
				}
				break;
			case 2:
				printf("Chemin absolu du répertoire : ");
				scanf("%s", g_dataFilePath);
				printf("\nNuméro du Repeat Cycle : ");
				scanf("%d", &RC);
				printf("Nom du canal : ");
				scanf("%s", nameChannel);
				printf("\nRésolution de l'image : ");
				scanf("%d", &resolution);
				printf("\n");
				if (TestParam(nameChannel,resolution) == 1){
					printf("... décodage en cours ...\n");
					if (strcmp(nameChannel,"all") == 0)
						AssembleAllImages(RC,0,70,resolution);
					else
						AssembleImages(RC,0,70,nameChannel,resolution);
				}
				break;
			case 3:
				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				if (TestFile(g_dataFile) != 1)
					break;
				printf("\nNuméro du Repeat Cycle : ");
				scanf("%d", &RC);
				printf("\nNuméro du swath : ");
				scanf("%d", &swath);
				printf("\n");
				printf("... décodage en cours ...\n");
				WriteScanAngles(swath,RC);
				break;
			case 4:
				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				if (TestFile(g_dataFile) != 1)
					break;
				printf("Numéro du swath : ");
				scanf("%d", &swath);
				printf("\n");
				printf("... décodage en cours ...\n");
				GetRCStartTime(parameterMOF);
				break;
			default:
				printf("Vous n'avez pas rentré un nombre correct.\n");
				break;
			}

			break;

			case 2:
				strcpy(tmp2,g_dataModelPath);
				strcat(tmp2,"MODEL_LI.eas");
				strcpy(tmp1,g_dataFilePath);
				strcat(tmp1,"LI_files/");
				g_dataFilePath = tmp1;
				g_dataModel = tmp2;

				status = select_DDR_i(g_dataModel);
				Check(status,__LINE__-1,__FILE__);

				int choixLI = 0;
				printf("==========================    LI    ==========================\n\n");
				printf("[0] Full Disk [1] OC 1  [2] OC 2  [3] OC 3  [4] OC 4  [5] Time\n\n");
				scanf("%d", &choixLI);
				if ((choixLI != 0) && (choixLI != 1) && (choixLI != 2) && (choixLI != 3) && (choixLI != 4)&& (choixLI != 5)){
					printf("Vous n'avez pas rentré un nombre correct.\n");
					break;
				}
				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				printf("\n");
				if (TestFile(g_dataFile) != 1)
					break;
				printf("\nNuméro du Repeat Cycle : ");
				scanf("%d", &RC);
				printf("\n");
				printf("... décodage en cours ...\n");
				if (choixLI==5){
					status =  TimeLI(parameterMOF);
					Check(status,__LINE__-1,__FILE__);
				}else{
					GetImageLI(RC, choixLI);
				}
				break;
			case 3:
				strcpy(tmp2,g_dataModelPath);
				strcat(tmp2,"MODEL_header.eas");
				g_dataModel = tmp2;
				status = select_DDR_i(g_dataModel);
				Check(status,__LINE__-1,__FILE__);

				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				if (TestFile(g_dataFile) != 1)
					break;
				printf("\nbloc : ");
				scanf("%d", &numeroBloc);
				printf("\n");

				printf("... décodage en cours ...\n");
				status = Header(parameterMOF, numeroBloc);
				Check(status,__LINE__-1,__FILE__);
				break;
			case 4:

				strcpy(tmp2,g_dataModelPath);
				strcat(tmp2,"MODEL_plateform.eas");
				g_dataModel = tmp2;
				status = select_DDR_i(g_dataModel);
				Check(status,__LINE__-1,__FILE__);

				printf("Chemin absolu du fichier : ");
				scanf("%s", g_dataFile);
				if (TestFile(g_dataFile) != 1)
					break;
				printf("... décodage en cours ...\n\n");

				Platform(parameterMOF);

				quitter = 0;
				break;
			case 5:
				quitter = 0;

				break;
			default:
				printf("Vous n'avez pas rentré un nombre correct.\n");
				break;
		}

		time(&arrivee);
		int secondes = ((int) difftime(arrivee, depart)) %60;
		int minutes = ((int) difftime(arrivee, depart)) /60;
		printf("\ntemps écoulé : %d m %d s.\n", minutes,secondes);


		status = stop_i();  //End of services EAST
		Check(status,__LINE__-1,__FILE__);

	}



	return EXIT_SUCCESS;
}






