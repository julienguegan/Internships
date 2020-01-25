/*
 * paquetReading.c
 *
 * This methods are created in purpose to get specifics data (like Scan angle, values of pixel...)
 * of a binary file by using East library
 *
 *  Created on: Apr 6, 2018
 *      Author: datamgr
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "interpreter.h"
#include "paquetReading.h"
#include "data.h"
#include "image.h"


const corr CHANNEL_RESOLUTION[17] = {{"VIS_1",224},{"VIS_2",224},{"VIS_3",224},{"VIS_4",224},{"VIS_5",448},{"NIR_1",224},{"NIR_2",224},{"NIR_3",448},
		{"IR1_1",224},{"IR1_2",224},{"IR1_3",112},{"IR1_4",112},{"IR2_1",112},{"IR2_2",112},{"IR3_1",224},{"IR3_2",112},{"IR3_3",112}} ;

/* - Get the Scan angles
 *
 * - input : a *SCAaxis = an array of SCAaxis which each element contain 50 flag and angles North/South
 *			 a *SCAaxis = an array of SCAaxis which each element 50 flag and angles East/West
 *			 an **int = a 2D array
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetScanAngles(SCAaxis *a_NSaxis,SCAaxis *a_EWaxis,int **a_fineTime,int *nbrSCA)
{
	char cible[] = "PDF.SDF.SCA.SCAN_OPERATION.BLOCK";
	int cptSCA = 0;
	int subtype,status;
	OCTET *buffer;
	int longueur;
	int i;
	int paquetNumber=0;
	while (1) { 	// Loop on the number of Block = occurence of our data model
		paquetNumber++;
		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status!=1)	{//stop there is an error
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}

		status = prepare_accesses_for_current_block_i();
		if (status!=1)	{//stop there is an error
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");		//get the Subtype

		if (subtype == 31) { 	// if subt = 31, it is a Scan Operation packet

			for (i=0; i< 50 ; i++) { 	//loop on our 50 angles

				a_fineTime[cptSCA][i] = GetIntInArray(cible, "FINE_TIME_TAG",i);		//get the fine time tag

				a_NSaxis[cptSCA].angle[i] = GetIntInArray(cible, "NS_AXIS.ANGLE_NS",i);		// get the NS scan angle
				a_NSaxis[cptSCA].angle[i] = a_NSaxis[cptSCA].angle[i]*360/pow(2,25);

				a_EWaxis[cptSCA].angle[i] = GetIntInArray(cible, "EW_AXIS.ANGLE_EW",i);	// get the EW scan angle
				a_EWaxis[cptSCA].angle[i] = a_EWaxis[cptSCA].angle[i]*360/pow(2,25);

			}	//end for

/*			if (cptSCA ==0) {
				i = 0;
				printf("\nflag   int      hexa      logical \n");
				printf("%d   %d ",GetIntInArray(cible, "NS_AXIS.FLAG_NS",i),GetIntInArray(cible, "NS_AXIS.ANGLE_NS",i));
				printf("   %x    ",(unsigned int) GetIntInArray(cible, "NS_AXIS.ANGLE_NS",i));
				char chaine[300];char temp[50];int k;
				strcpy(chaine,cible);
				sprintf(temp,"(%d).",i+1);
				strcat(chaine,temp);
				strcat(chaine,"NS_AXIS.ANGLE_NS");
				status = get_data_entity_i(chaine,&buffer,&longueur,LOGICAL);
				CheckEastStatus(status);
				for (k=0;k<sizeof(buffer);k++)
					printf("%x",buffer[k]);
				printf("\n");
				free(buffer);
			}*/

			cptSCA++;

		}	//end if

	}	//end while

	*nbrSCA = cptSCA;
	return 1;
}

/* - Get all pixel of a channel which can appear in 1, 2 or 3 types of pixel data raw frame packets
 * 	 according to the chosen channel
 *
 * - input : an **int = a 2D array where the value of the pixel are stored (dimensions : nbrPixelByColumn*nbrPixelByRow)
 * 			 a string = name of the channel we want to get back (ex : VIS_5, NIR_3)
 * 			 an *int = number of Pixel by Row => width of the image to be updated, needed for creation of the image
 * 			 an int = indicate the resolution of the chosen channel, example : 3 means pick one pixel out of three (in the 2 directions)
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetChannel(int **a_channel,char *a_nameChannel, int *nbrPixelByRow, int n)
{
	char cible[] = "PDF.SDF.VCU.PIXEL_DATA";
	int status, subtype, i = 0, j = 0, k;
	int size ;
	for (k=0;k<17;k++){
		if (strcmp(CHANNEL_RESOLUTION[k].channel,a_nameChannel) == 0)
			size = CHANNEL_RESOLUTION[k].size;
	}

	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status!=1)	{//stop there is an error
			printf("Erreur East %d  : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status!=1)	{//stop there is an error
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");  //get the subtype

		if (subtype == 3) {

			if ( j % n == 0 ) {

				char tmp[200];
				strcpy(tmp,cible);
				strcat(tmp,".RAW_FRAME_4."); //create the full path to data
				strcat(tmp,a_nameChannel);

				for (i=0; i< size/n ; i++) 	 //loop on our bytes = one column for our image
					a_channel[i][j/n] = GetIntArray(tmp,n*i);

			}
			j++;

		} else if ((subtype == 2) && ( (size == 224) || (size == 448) ) ){  //raw frame 4

			if ( j % n == 0 ) {

				char tmp[200];
				strcpy(tmp,cible);
				strcat(tmp,".RAW_FRAME_2.");
				strcat(tmp,a_nameChannel);

				for (i=0; i< size/n ; i++)
					a_channel[i][j/n] = GetIntArray(tmp,n*i);

			}
			j++;

		} else if ( (subtype == 1) && (size == 448) ){

			if ( j % n == 0 ) {

				char tmp[200];
				strcpy(tmp,cible);
				strcat(tmp,".RAW_FRAME_1.");
				strcat(tmp,a_nameChannel);

				for (i=0; i< size/n ; i++)
					a_channel[i][j/n] = GetIntArray(tmp,n*i);

			}

			j++;

		}

	}	//end while

	*nbrPixelByRow = j;

	return 1;

}

/* - Get all pixel of all channel in 3D array
 *
 * - input : an **int = a 2D array where the value of the pixel are stored (dimensions : nbrPixelByColumn*nbrPixelByRow)
 * 			 an int = indicate the resolution of the chosen channel, example : 3 means pick one pixel out of three (in the 2 directions)
 * 			 an int = number of Pixel by Row for raw frame 1&3
 * 			 an int = number of Pixel by Row for raw frame 2
 * 			 an int = number of Pixel by Row for raw frame 4
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetAllChannel(int ***a_channel, int n, int *taille1, int *taille2, int *taille3)
{
	char cible[] = "PDF.SDF.VCU.PIXEL_DATA";
	int status, subtype, k, i, size, j1 = -1, j2 = -1, j3 = -1;
	char *nameChannel;
	char tmp[200];
	char nameRawFrame[50] = ".RAW_FRAME_x.";

	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status == -7)	//stop if it is the end of the file
			break;
		if (status!=1)	{//stop there is an error
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status!=1)	{
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");  //get the subtype

		if ((subtype == 1) || (subtype == 2) || (subtype == 3)) {

			if (subtype == 1) {
				nameRawFrame[11] = '1';
				j1++;

			} else if (subtype == 2) {
				nameRawFrame[11] = '2';
				j2++;
				j1++;
			} else if (subtype == 3) {
				nameRawFrame[11] = '4';
				j3++;
				j2++;
				j1++;
			}

			for (k=0;k<17;k++){

				if (CHANNEL_RESOLUTION[k].size >= (int) (448/pow(2,subtype-1))) { //les channel qui sont dans ce subtype

					size = CHANNEL_RESOLUTION[k].size;
					nameChannel = CHANNEL_RESOLUTION[k].channel;
					strcpy(tmp,cible);
					strcat(tmp,nameRawFrame);
					strcat(tmp,nameChannel);

					if (size == 448) {
						if (j1 % n == 0) {
							for (i=0; i< size/n ; i++)
								a_channel[k][i][j1/n] = GetIntArray(tmp,n*i);
						}
					} else if (size == 224) {
						if (j2 % n == 0) {
							for (i=0; i< size/n ; i++)
								a_channel[k][i][j2/n] = GetIntArray(tmp,n*i);
						}
					} else if (size == 112) {
						if (j3 % n == 0) {
							for (i=0; i< size/n ; i++)
								a_channel[k][i][j3/n] = GetIntArray(tmp,n*i);
						}
					}
				}
			}

		}
	}
	*taille1=j1;
	*taille2=j2;
	*taille3=j3;
	return 1;
}

/*	- Save image of an OC by getting all pixel value
 *
 * - input : an int (1, 2, 3 or 4) = correspond to the Optical Channel (=PID) wanted
 * 			 a string = name of the image
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetOC(int OC, char nameImageOC[100])
{
	int status, subtype, i, j, k, m, n, p, q, w, imageWINDOW[10][13], win_counter, ASIC_ID, N_packet;
	int **imageOC;
	imageOC = malloc(1000 * sizeof(*imageOC));
	if (imageOC == NULL){
		printf("erreur d'allocation : ligne %d du fichier %s \n",__LINE__,__FILE__);
		return -1;
	}
	for(i=0 ; i < 1000 ; i++) {
		imageOC[i] = malloc(1170 * sizeof(**imageOC));
		if (imageOC[i] == NULL) {
			printf("erreur d'allocation : ligne %d du fichier %s \n",__LINE__,__FILE__);
			return -1;
		}
	}

	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		if (GetInt("PH.PID.APID.PID")-63 == OC){ //the PID = OC # (changed 64/65/66/67 in 1/2/3/4)

			subtype = GetInt("PDF.DFH.SUBT");  //get the subtype

			if (subtype == 63) {  // 63 =  img_data packet

				N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");

				for (k=0; k<N_packet; k++){	// number of window in the actual packet

					ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;

					win_counter = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.WIN_COUNTER", k);	//2250 win_counter

					for (w=0; w<130; w++){ // 10*13 window reference

						m = w / 13; // row n° of a WINDOW : de 0 to 9
						n = w % 13; // col n° of a WINDOW : de 0 to 12

						imageWINDOW[m][n] = GetIntIn2Array("PDF.SDF.IMG_DATA.REPEATED", "BKG_PIX", k, w, "BACKGROUND_VALUE_PIXEL");

						q = n + (win_counter * 13) % 585;	//  col n° of a ASIC =  0 to 585
						p = m + (win_counter * 13) / 585 * 10;	// row n° of a ASIC = 0 to 500

						j = p + m; // row n° of a OC = 0 to 1000
						i = q + n; // col n° of a OC = 0 to 1170

						if (ASIC_ID == 1) {
							i = q;
							j = p;
						} else if (ASIC_ID == 2) {
							i = 585 + q;
							j = p;
						} else if (ASIC_ID == 3)  {
							i = q;
							j = 999 - p;
						} else if (ASIC_ID == 4) {
							i = 585 + q;
							j = 999 - p;
						}

						imageOC[j][i] = imageWINDOW[m][n];

					}

				}

			}

		}

	}

	CreateImage(imageOC,nameImageOC,1170,1000,"noflip");

	for(i=0 ; i < 1000 ; i++) {
		free(imageOC[i]);
	}
	free(imageOC);

	return 1;

}

/*	- Save full image of the LI (= 4 OC images ) by getting all pixel value
 *
 * - input : a string : name of the image
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 *
 */
int GetFullOC(char nameImage[100])
{
	int status, subtype, i, j, k, m, n, p, q, w, x, y, win_counter, ASIC_ID, N_packet, OC;
	int imageOC[1000][1170] , imageWINDOW[10][13] ;
	int **imageFull;
	imageFull = malloc(2000 * sizeof(*imageFull));
	for(i=0 ; i < 2000 ; i++)
		imageFull[i] = malloc(2340 * sizeof(**imageFull));

	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		OC = GetInt("PH.PID.APID.PID")-63; //the PID = OC # (changed 64/65/66/67 in 1/2/3/4)

		subtype = GetInt("PDF.DFH.SUBT");  //get the subtype

		if (subtype == 63) {  // 63 =  img_data packet

			N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");

			for (k=0; k<N_packet; k++){	// number of window in the actual packet

				ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;

				win_counter = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.WIN_COUNTER", k);	//2250 win_counter

				for (w=0; w<130; w++){ // 10*13 window reference

					m = w / 13; // row n° of a WINDOW : de 0 to 9
					n = w % 13; // col n° of a WINDOW : de 0 to 12

					imageWINDOW[m][n] = GetIntIn2Array("PDF.SDF.IMG_DATA.REPEATED", "BKG_PIX", k, w, "BACKGROUND_VALUE_PIXEL");

					q = n + (win_counter * 13) % 585;		//  col n° of a ASIC =  0 to 585
					p = m + (win_counter * 13) / 585 * 10;	// row n° of a ASIC = 0 to 500

					if (ASIC_ID == 1) {
						i = q;					    // col n° of a OC = 0 to 1170
						j = p;  					// row n° of a OC = 0 to 1000
					} else if (ASIC_ID == 2) {
						i = 585 + q;
						j = p;
					} else if (ASIC_ID == 3) L {
						i = q;
						j = 999 - p;
					} else if (ASIC_ID == 4) {
						i = 585 + q;
						j = 999 - p;
					}

					imageOC[j][i] = imageWINDOW[m][n];

					if (OC == 1) {
						x = 1000 + j;
						y = 1170 + i;
					} else if (OC == 2) {
						x = 1000 + j;
						y = i;
					} else if (OC == 3)  {
						x = 999 - j;
						y = 1169 - i;
					} else if (OC == 4) {
						x = 999 - j;
						y = 2339 - i;
					}

					imageFull[x][y] = imageOC[j][i] ;

				}

			}

		}

	}

	CreateImage(imageFull, nameImage, 2340, 2000, "noflip");


	for(i=0 ; i < 2000 ; i++)
		free(imageFull[i]);
	free(imageFull);

	return 1;
}

int TimeLI(parameterMOF parameterMOF) {

	int i, j, status, OC_ID, c1, c2, c3, c4, f1, f2, f3, k1=0, k2=0, k3=0, k4=0, k=0;
	double OBT;
	double **timeOC;
	timeOC = malloc(100000 * sizeof(*timeOC));
	for(i=0 ; i < 100000 ; i++)
		timeOC[i] = malloc(100000 * sizeof(**timeOC));
	double *timeallOC = malloc(100000 * sizeof(*timeallOC));

	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		c1 = GetInt("PDF.DFH.TM_TIME.C1");
		c2 = GetInt("PDF.DFH.TM_TIME.C2");
		c3 = GetInt("PDF.DFH.TM_TIME.C3");
		c4 = GetInt("PDF.DFH.TM_TIME.C4");
		f1 = GetInt("PDF.DFH.TM_TIME.F1");
		f2 = GetInt("PDF.DFH.TM_TIME.F2");
		f3 = GetInt("PDF.DFH.TM_TIME.F3");
		OBT = c1*pow(256,3)+c2*pow(256,2)+c3*pow(256,1)+c4*pow(256,0)+f1*pow(256,-1)+f2*pow(256,-2)+f3*pow(256,-3);

		OC_ID = GetInt("PH.PID.APID.PID")-63; //the PID = OC # (changed 64/65/66/67 in 1/2/3/4)
		timeallOC[k] = OBT;
		k++;
		if (OC_ID == 1){
			timeOC[0][k1] = OBT;
			k1 = k1+1;
		}else if (OC_ID == 2){
			timeOC[1][k2] = OBT;
			k2 = k2+1;
		}else if (OC_ID == 3){
			timeOC[2][k3] = OBT;
			k3 = k3+1;
		}else if (OC_ID == 4){
			timeOC[3][k4] = OBT;
			k4 = k4+1;
		}

	}

	FILE* fichierEcrire = NULL;
	fichierEcrire = fopen("timeOC.csv","w");
	if (fichierEcrire != NULL) {
		printf("écriture en cours ... \n");
		fprintf(fichierEcrire,"OC1 OC2 OC3 OC4 \n");
		for (j=0;j<k1;j++)
			fprintf(fichierEcrire,"%lf %lf %lf %lf \n",timeOC[0][j],timeOC[1][j],timeOC[2][j],timeOC[3][j]);
		fclose(fichierEcrire);
		printf("fichier écrit !\n");
	}

		fichierEcrire = fopen("timeallOC.csv","w");
		if (fichierEcrire != NULL) {
			printf("écriture en cours ... \n");
			for (j=0;j<k;j++)
				fprintf(fichierEcrire,"%lf \n",timeallOC[j]);
			fclose(fichierEcrire);
			printf("fichier écrit !\n");
		}


	for(i=0 ; i < 100000 ; i++)
		free(timeOC[i]);
	free(timeOC);
	free(timeallOC);

	return 1;
}

int TimeLIasic(parameterMOF parameterMOF) {

	int i, j, k, status, OC_ID, ASIC_ID, N_packet, c1, c2, c3, c4, f1, f2, f3, k1=0, k2=0, k3=0, k4=0;
	double OBT,UTC;
	double ***timeOCasic;
	timeOCasic = malloc(4 * sizeof(*timeOCasic));
	for(i=0 ; i < 4 ; i++)
		timeOCasic[i] = malloc(4 * sizeof(**timeOCasic));
	for(i=0 ; i < 4 ; i++){
		for(j=0 ; j < 4 ; j++)
			timeOCasic[i][j] = malloc(100000 * sizeof(***timeOCasic));
	}


	while (1) {		 // Loop on the number of occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		c1 = GetInt("PDF.DFH.TM_TIME.C1");
		c2 = GetInt("PDF.DFH.TM_TIME.C2");
		c3 = GetInt("PDF.DFH.TM_TIME.C3");
		c4 = GetInt("PDF.DFH.TM_TIME.C4");
		f1 = GetInt("PDF.DFH.TM_TIME.F1");
		f2 = GetInt("PDF.DFH.TM_TIME.F2");
		f3 = GetInt("PDF.DFH.TM_TIME.F3");
		OBT = c1*pow(256,3)+c2*pow(256,2)+c3*pow(256,1)+c4*pow(256,0)+f1*pow(256,-1)+f2*pow(256,-2)+f3*pow(256,-3);
		UTC = parameterMOF.UTC0 + parameterMOF.gradient * (OBT - parameterMOF.OBT0) + parameterMOF.offset;

		OC_ID = GetInt("PH.PID.APID.PID")-63; //the PID = OC # (changed 64/65/66/67 in 1/2/3/4)

		switch (OC_ID)  {
		case 1:
			if (GetInt("PDF.DFH.SUBT") == 63) {  // 63 =  img_data packet
				N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");
				for (k=0; k<N_packet; k++){  // number of window in the actual packet
					ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;
					timeOCasic[OC_ID-1][ASIC_ID-1][k1] = UTC;
					k1 = k1+1;
				}
			}
			break;
		case 2:

			if (GetInt("PDF.DFH.SUBT") == 63) {  // 63 =  img_data packet
				N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");
				for (k=0; k<N_packet; k++){  // number of window in the actual packet
					ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;
					timeOCasic[OC_ID-1][ASIC_ID-1][k2] = UTC;
					k2 = k2+1;
				}
			}
			break;
		case 3:

			if (GetInt("PDF.DFH.SUBT") == 63) {  // 63 =  img_data packet
				N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");
				for (k=0; k<N_packet; k++){  // number of window in the actual packet
					ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;
					timeOCasic[OC_ID-1][ASIC_ID-1][k3] = UTC;
					k3 = k3+1;
				}
			}
			break;
		case 4:
			if (GetInt("PDF.DFH.SUBT") == 63) {  // 63 =  img_data packet
				N_packet = GetInt("PDF.SDF.IMG_DATA.N_PACKETS");
				for (k=0; k<N_packet; k++){  // number of window in the actual packet
					ASIC_ID = GetIntInArray("PDF.SDF.IMG_DATA.REPEATED", "WINDOWS_IDENTIFICATION.ASIC_ID", k)+1;
					timeOCasic[OC_ID-1][ASIC_ID-1][k4] = UTC;
					k4 = k4+1;
				}
			}
			break;

		default:
			printf("problemo.\n");
			break;
		}

	}

	FILE* fichierEcrire = NULL;
	fichierEcrire = fopen("timeOCasic.csv","w");
	if (fichierEcrire != NULL) {
		printf("écriture en cours ... \n");
		fprintf(fichierEcrire,"OC1... OC2... OC3... OC4... \n");
		fprintf(fichierEcrire,"ASIC1 ASIC2 ASIC3 ASIC4 \n");
		for (j=0;j<k1;j++){
			fprintf(fichierEcrire,"%lf %lf %lf %lf ",timeOCasic[0][0][j],timeOCasic[0][1][j],timeOCasic[0][2][j],timeOCasic[0][3][j]);
			fprintf(fichierEcrire,"%lf %lf %lf %lf ",timeOCasic[1][0][j],timeOCasic[1][1][j],timeOCasic[1][2][j],timeOCasic[1][3][j]);
			fprintf(fichierEcrire,"%lf %lf %lf %lf ",timeOCasic[2][0][j],timeOCasic[2][1][j],timeOCasic[2][2][j],timeOCasic[2][3][j]);
			fprintf(fichierEcrire,"%lf %lf %lf %lf \n",timeOCasic[3][0][j],timeOCasic[3][1][j],timeOCasic[3][2][j],timeOCasic[3][3][j]);
		}
		fclose(fichierEcrire);
		printf("fichier écrit !\n");
	}

	for(k=0 ; k < 4 ; k++){
		for(i=0 ; i < 4 ; i++)
			free(timeOCasic[k][i]);
		free(timeOCasic[k]);
	}
	free(timeOCasic);

	return 1;
}

/* - Print the date of the RC start time in human readable format
 *
 * - input : char *date = name of the image
 *  		 long double = UTC0
 *			 int = gradient
 *			 long double = OBT0
 *			 int  = offset
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetRCStartTime(parameterMOF parameterMOF ){

	int subtype, status, c1, c2, c3, c4, f1, f2, f3;
	double OBT, UTC;
	int test = 0;
	while (1) { 	// Loop on the number of Block = occurence of our data model

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file
		if (status==-7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");  //get the subtype

		if (subtype == 41) {

			c1 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.C1");
			c2 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.C2");
			c3 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.C3");
			c4 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.C4");
			f1 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.F1");
			f2 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.F2");
			f3 = GetInt("PDF.SDF.ICU.STATIC_INFO.REPEAT_CYCLE_START_TIME.F3");
			test = 1;
			break;
		}
	}
	if (test == 1) {
		OBT = c1*pow(256,3)+c2*pow(256,2)+c3*pow(256,1)+c4*pow(256,0)+f1*pow(256,-1)+f2*pow(256,-2)+f3*pow(256,-3);
		UTC = parameterMOF.UTC0 + parameterMOF.gradient * (OBT - parameterMOF.OBT0) + parameterMOF.offset;
		printf("OBT = %.6f ( %x %x %x %x %x %x %x ) \nUTC = %.6f \n",OBT,c1,c2,c3,c4,f1,f2,f3, UTC);
		time_t utc2;
		utc2 = (time_t)UTC;
		printf("RC start time : %s\n",asctime(gmtime(&utc2)));
	} else
		printf("Pas de paquet ICU dans ce fichier.\n");

	return 1;
}

/* - Print useful informations of the header of a packet
 *
 * - input : char *date = name of the image
 *  		 long double = UTC0
 *			 int = gradient
 *			 long double = OBT0
 *			 int  = offset
 *
 * - output : print on the console
 * 			  int = 1 if no error , -1 if error
 */
int Header(parameterMOF parameterMOF, int numeroBloc){

	int status, apid, bloc = 0, c1, c2, c3, c4, f1, f2, f3, subtype, ver, cat, dfhf, grpf, seqcnt, plen;
	double OBT, UTC;
	char *date;

	while (1) {

		bloc++;

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file

		if (status == -7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}

		if (bloc == numeroBloc) {

			ver = GetInt("PH.PID.VER");
			cat = GetInt("PH.PID.CAT");
			dfhf = GetInt("PH.PID.DFHF");
			apid = GetInt("PH.PID.APID");
			grpf = GetInt("PH.PSC.GRPF");
			seqcnt = GetInt("PH.PSC.SEQCNT");
			plen = GetInt("PH.PLEN");

			subtype = GetInt("PDF.DFH.SUBT");

			c1 = GetInt("PDF.DFH.TIME.C1");
			c2 = GetInt("PDF.DFH.TIME.C2");
			c3 = GetInt("PDF.DFH.TIME.C3");
			c4 = GetInt("PDF.DFH.TIME.C4");
			f1 = GetInt("PDF.DFH.TIME.F1");
			f2 = GetInt("PDF.DFH.TIME.F2");
			f3 = GetInt("PDF.DFH.TIME.F3");

			break;
		}

	}
	OBT = c1*pow(256,3)+c2*pow(256,2)+c3*pow(256,1)+c4*pow(256,0)+f1*pow(256,-1)+f2*pow(256,-2)+f3*pow(256,-3);
	UTC = parameterMOF.UTC0 + parameterMOF.gradient * (OBT - parameterMOF.OBT0) + parameterMOF.offset;
	time_t utc2;
	utc2 = (time_t)UTC;
	date = asctime(gmtime(&utc2));//ctime(&utc2);


	printf("bloc : %d  \n",bloc);
	printf("version : %d  \n",ver);
	printf("type: %d  \n",cat);
	printf("data field header flag : %d  \n",dfhf);
	printf("APID : %x  \n",apid);
	printf("grouping flags : %d  \n",grpf);
	printf("source sequence count : %d  \n",seqcnt);
	printf("packet length: %d  \n",plen);
	printf("subtype: %d  \n",subtype);
	printf("Absolute time : %s  \n",date);
	printf("OBT = %.6f ( %x %x %x %x %x %x %x ) \nUTC = %.6f \n",OBT,c1,c2,c3,c4,f1,f2,f3,UTC);

	return 1;
}

/* Create a file where we have the value of quaternion in a file and date (secondes and calendar fomat) in a csv file
 *
 * input : - MOF coefficient
 */
int Platform(parameterMOF parameterMOF){

	int status, bloc = 0, c1, c2, c3, c4, f1, f2, f3, subtype, k = 0;
	double OBT, q1, q2, q3, q4;
	int i;
	time_t utc2;

	double **quaternion = malloc(100000 * sizeof(*quaternion));
	for(i=0 ; i < 100000 ; i++)
		quaternion[i] = malloc(4 * sizeof(**quaternion));
	double *UTC = malloc(100000 * sizeof(UTC));

	char **date = malloc(10000 * sizeof(*date));
	for(i=0 ; i < 10000 ; i++)
		date[i] = malloc(100 * sizeof(**date));

	while (1) {

		bloc++;

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file

		if (status == -7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}

		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");

		if (subtype == 221) {

			c1 = GetInt("PDF.SDF.FPM.TIME.C1");
			c2 = GetInt("PDF.SDF.FPM.TIME.C2");
			c3 = GetInt("PDF.SDF.FPM.TIME.C3");
			c4 = GetInt("PDF.SDF.FPM.TIME.C4");
			f1 = GetInt("PDF.SDF.FPM.TIME.F1");
			f2 = GetInt("PDF.SDF.FPM.TIME.F2");
			f3 = GetInt("PDF.SDF.FPM.TIME.F3");

			q1 = GetDouble("PDF.SDF.FPM.QI_1");
			q2 = GetDouble("PDF.SDF.FPM.QI_2");
			q3 = GetDouble("PDF.SDF.FPM.QI_3");
			q4 = GetDouble("PDF.SDF.FPM.QI_4");

			quaternion[k][0] = q1;
			quaternion[k][1] = q2;
			quaternion[k][2] = q3;
			quaternion[k][3] = q4;

			OBT = c1*pow(256,3)+c2*pow(256,2)+c3*pow(256,1)+c4*pow(256,0)+f1*pow(256,-1)+f2*pow(256,-2)+f3*pow(256,-3);
			UTC[k] = (double) parameterMOF.UTC0 + parameterMOF.gradient * (OBT - parameterMOF.OBT0) + parameterMOF.offset;

			utc2 = (time_t) UTC[k];
			strcpy(date[k],asctime(gmtime(&utc2)));

			k++;

		}

	}

	if (k == 0)
		printf("pas de décodage\n");

	FILE* fichierEcrire = NULL;
	fichierEcrire = fopen("quaternion.csv","w");
	if (fichierEcrire != NULL) {
		printf("écriture en cours ... \n");
		fprintf(fichierEcrire,"date secondes Q1 Q2 Q3 Q4\n");
		for (i=0;i<k;i++) { // lines
			fprintf(fichierEcrire,"%s %lf %lf %lf %lf %lf\n",date[i],UTC[i],quaternion[i][0],quaternion[i][1],quaternion[i][2],quaternion[i][3]);
		}
		fclose(fichierEcrire);
		printf("fichier écrit !\n");
	}

	return 1;
}

/* - Count the number of the differents raw frame packets in our file by fully browsing it
 *
 * - input : an *int = the number of raw frame 1&3 packets
 * 			 an *int = the number of raw frame 2 packets
 * 			 an *int = the number of raw frame 4 packets
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetNumberOfFrame(int *a_nbrFrame13,int *a_nbrFrame2,int *a_nbrFrame4){

	int status = 0, subtype = 0;

	while (1) {


		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file

		if (status == -7)	//stop if it is the end of the file
			break;
		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}
		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");

		if (subtype == 1)
			*a_nbrFrame13 = *a_nbrFrame13+1;
		else if (subtype == 2)
			*a_nbrFrame2 = *a_nbrFrame2+1;
		else if (subtype == 3)
			*a_nbrFrame4 = *a_nbrFrame4+1;



	}

	return 1;
}

/* - Count the number of the SCA packets in our file by fully browsing it
 *
 * - input : an *int = the number of SCA packets
 *
 * - output : update thanks to pointer
 * 			  int = 1 if no error , -1 if error
 */
int GetNumberOfSCA(int *a_nbrSCA){

	int status = 0, subtype = 0;

	while (1) {

		status = load_next_data_block_from_file_i(g_dataFile);	//read the data file

		if (status == -7)	//stop if it is the end of the file
			break;

		if (status != 1) {
			printf("Erreur East %d (fichier de données d'entrée) : ligne %d du fichier %s \n",status,__LINE__-4,__FILE__);
			return -1;
		}

		status = prepare_accesses_for_current_block_i();
		if (status != 1) {
			printf("Erreur East %d : ligne %d du fichier %s \n",status,__LINE__-2,__FILE__);
			return -1;
		}

		subtype = GetInt("PDF.DFH.SUBT");

		if (subtype == 31)
			*a_nbrSCA = *a_nbrSCA+1;

	}
	return 1;
}


