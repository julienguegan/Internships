/*
 * createImage.c
 *
 *  Created on: Apr 16, 2018
 *      Author: datamgr
 */

/*
 * 01_09_02_bmp.c
 *
 *  Created on: Apr 16, 2018
 *      Author: datamgr
 */


#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "image.h"

/* Initialize an image of dimension width*heigth
 *
 * input : an int = the width
 * 		   an int = the heigth
 *
 * output : an Image structure
 */
Image* NouvelleImage(int w, int h)
{
	Image* I = malloc(sizeof(Image));
	I->w = w;
	I->h = h;
	I->dat = calloc(1,w*h*sizeof(Pixel*));
	return I;
}

/* Deallocate an image
 *
 */
void DeleteImage(Image* I)
{
	if (I)
	{
		free(I->dat);
		free(I);
	}
}

/* Set a pixel p in an image I at the position (i,j)
 *
 * intput : an Image* = the image structure where we want to set our pixel
 * 			an int = the position along the horizontal axis
 * 			an int = the position along the vertical axis
 * 			a Pixel = the pixel to set
 */
void SetPixel(Image* I, int i, int j, Pixel p)
{
	assert(I && i>=0 && i<I->w && j>=0 && j<I->h);
	I->dat[I->w*j+i] = p;
}


/* Get a pixel p in an image I at the position (i,j)
 *
 * output : a Pixel
 */
Pixel GetPixel(Image* I, int i, int j)
{
	assert(I && i>=0 && i<I->w && j>=0 && j<I->h);
	return I->dat[I->w*j+i];
}

// -------------------------------------------

#pragma pack(1)  // desative l'alignement mémoire
typedef int int32;
typedef short int16;

struct BMPImHead
{
	int32 size_imhead;
	int32 width;
	int32 height;
	int16 nbplans; // toujours 1
	int16 bpp;
	int32 compression;
	int32 sizeim;
	int32 hres;
	int32 vres;
	int32 cpalette;
	int32 cIpalette;
};

struct BMPHead
{
	char signature[2];
	int32 taille;
	int32 rsv;
	int32 offsetim;
	struct BMPImHead imhead;
};

Image* Charger(const char* fichier)
{
	struct BMPHead head;
	Image* I;
	int i,j,pitch;
	unsigned char bgrpix[3];
	char corrpitch[4] = {0,3,2,1};
	Pixel p;
	FILE* F = fopen(fichier,"rb");
	if (!F)
		return NULL;
	fread(&head,sizeof(struct BMPHead),1,F);
	if (head.signature[0] != 'B'  ||  head.signature[1] != 'M')
		return NULL;  // mauvaise signature, ou BMP non supporté.
	if (head.imhead.bpp != 24)
		return NULL;  // supporte que le 24 bits pour l'instant
	if (head.imhead.compression != 0)
		return NULL;  // rarrissime, je ne sais même pas si un logiciel écrit/lit des BMP compressés.
	if (head.imhead.cpalette !=0  ||  head.imhead.cIpalette !=0 )
		return NULL; // pas de palette supportée, cependant, a bpp 24, il n'y a pas de palette.
	I = NouvelleImage(head.imhead.width,head.imhead.height);
	pitch = corrpitch[(3*head.imhead.width)%4];
	for(j=0;j<I->h;j++)
	{
		for(i=0;i<I->w;i++)
		{
			fread(&bgrpix,1,3,F);
			p.r = bgrpix[2];
			p.g = bgrpix[1];
			p.b = bgrpix[0];
			SetPixel(I,i,I->h-j-1,p);
		}
		fread(&bgrpix,1,pitch,F);
	}
	fclose(F);
	return I;
}

/* Save an image bitmap by filling the header and the pixel
 *
 * input : an image structure
 * 		   a string = name of the file
 *
 */
int SaveImage(Image* I, const char* fichier)
{
	struct BMPHead head;
	Pixel p;
	int i,j,tailledata,pitch;
	unsigned char bgrpix[3];
	char corrpitch[4] = {0,3,2,1};
	FILE* F = fopen(fichier,"wb");
	if (!F)
		return -1;
	memset(&head,0,sizeof(struct BMPHead));
	head.signature[0] = 'B';
	head.signature[1] = 'M';
	head.offsetim = sizeof(struct BMPHead); // je vais toujours écrire sur le même moule.
	head.imhead.size_imhead = sizeof(struct BMPImHead);
	head.imhead.width = I->w;
	head.imhead.height = I->h;
	head.imhead.nbplans = 1;
	head.imhead.bpp = 24;
	pitch = corrpitch[(3*head.imhead.width)%4];
	tailledata = 3*head.imhead.height*head.imhead.width + head.imhead.height*pitch;
	head.imhead.sizeim = tailledata;
	head.taille = head.offsetim + head.imhead.sizeim;
	fwrite(&head,sizeof(struct BMPHead),1,F);
	for(j=0;j<I->h;j++)
	{
		for(i=0;i<I->w;i++)
		{
			p = GetPixel(I,i,I->h-j-1);
			bgrpix[0] = p.b;
			bgrpix[1] = p.g;
			bgrpix[2] = p.r;
			fwrite(&bgrpix,1,3,F);
		}
		bgrpix[0] = bgrpix[1] = bgrpix[2] = 0;
		fwrite(&bgrpix,1,pitch,F);
	}
	fclose(F);
	return 0;
}

/* Add an image to a bigger one
 *
 * - input : Image* = image to add
 * 			 Image* = bigger image
 * 			 int = x-axis position, where to start to add the image
 * 			 int = y-axis position, where to start to add the image
 * 			 char* = indicate if we have to flip or not the image
 * 			 int = size of the Uturn to not add
 */
void AddImageNoUturn(Image* image, Image* imageAssemblee, int x, int y, char* flip, int Uturn){

	int i,j;
	if (strcmp(flip,"flip") == 0){

		for(i=Uturn;i<image->w;i++) { //start at a certain pixel to not read the U-turn
			for(j=0;j<(image->h);j++){
				Pixel p;
				p.r = GetPixel(image,i,j).r;
				p.g = GetPixel(image,i,j).g;
				p.b = GetPixel(image,i,j).b;

				SetPixel(imageAssemblee,i+x-Uturn,j+y,p);
			}
		}


	} else if (strcmp(flip,"noflip") == 0) {

		for(i=0;i<image->w-Uturn;i++) { //end at image.w-Uturn to not read the U-turn
			for(j=0;j<(image->h);j++){
				Pixel p;
				p.r = GetPixel(image,i,j).r;
				p.g = GetPixel(image,i,j).g;
				p.b = GetPixel(image,i,j).b;

				SetPixel(imageAssemblee,i+x,j+y,p);
			}
		}

	}

}

/* Save an image bitmap
 *
 * - input : int** = a 2D array containing the value of the pixel
 * 			 char* = name of the image
 * 			 int = width of the image
 * 			 int = height of the image
 * 			 char* = indicate if we have to flip or not the image
 */
void CreateImage(int** channel, char* nomImage, int nbrPixelByRow, int nbrPixelByColumn, char* flip){

	Image* monImageACreer = NouvelleImage(nbrPixelByRow,nbrPixelByColumn);//creation image : width*height = nbrVCU*nbrPixelColumn px
	int i, j, pixelValue;

	for(j=0;j<nbrPixelByRow;j++){	// along horizontal axis

		for(i=0;i<nbrPixelByColumn;i++){  // along vertical axis

			Pixel p;

			pixelValue = ((float)channel[i][j])*256/4096;	// = (2⁸/2¹²) convert a value of 12 bits in 8 bits

			p.r = pixelValue; // we want only 1 wavelength = no rgb = all at the same value
			p.g = pixelValue;
			p.b = pixelValue;

			if (strcmp(flip,"noflip") == 0)
				SetPixel(monImageACreer,j,i,p);
			else if (strcmp(flip,"flip horizontal") == 0)
				SetPixel(monImageACreer,nbrPixelByRow-j-1,i,p);
			else if (strcmp(flip,"flip vertical") == 0)
				SetPixel(monImageACreer,j,nbrPixelByColumn-i-1,p);
			else if (strcmp(flip,"flip horizontal vertical") == 0)
				SetPixel(monImageACreer,nbrPixelByRow-j-1,nbrPixelByColumn-i-1,p);
		}
	}
	SaveImage(monImageACreer,nomImage);
	printf("\"%s\" sauvegardée \n", nomImage);
	DeleteImage(monImageACreer);
}



