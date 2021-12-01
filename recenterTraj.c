#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct dump
{
	int sino, type, ix, iy, iz;
	float x, y, z, xs, ys, zs;
} DUMP;

DUMP findCOM (DUMP *atomEntries, int nAtoms)
{
	DUMP com;
	com.x = 0; com.y = 0; com.z = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		com.x += atomEntries[i].x;
		com.y += atomEntries[i].y;
		com.z += atomEntries[i].z;
	}

	com.x /= nAtoms; com.y /= nAtoms; com.z /= nAtoms;

	return com;
}

DUMP *recenterCoordinates (DUMP *atomEntries, int nAtoms, DUMP com)
{
	DUMP *recenteredEntries;
	recenteredEntries = (DUMP *) malloc (nAtoms * sizeof (DUMP));

	for (int i = 0; i < nAtoms; ++i)
	{
		recenteredEntries[i].x = atomEntries[i].x - com.x;
		recenteredEntries[i].y = atomEntries[i].y - com.y;
		recenteredEntries[i].z = atomEntries[i].z - com.z;
		recenteredEntries[i].sino = atomEntries[i].sino;
		recenteredEntries[i].type = atomEntries[i].type;
		recenteredEntries[i].ix = 0;
		recenteredEntries[i].iy = 0;
		recenteredEntries[i].iz = 0;
		// recenteredEntries[i].ix = atomEntries[i].ix;
		// recenteredEntries[i].iy = atomEntries[i].iy;
		// recenteredEntries[i].iz = atomEntries[i].iz;
	}

	return recenteredEntries;
}

void printCoordinates (DUMP *recenteredEntries, int nAtoms, FILE *output)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		fprintf(output, "%d %d %f %f %f %d %d %d\n", 
			recenteredEntries[i].sino, 
			recenteredEntries[i].type, 
			recenteredEntries[i].x, 
			recenteredEntries[i].y, 
			recenteredEntries[i].z, 
			recenteredEntries[i].ix, 
			recenteredEntries[i].iy, 
			recenteredEntries[i].iz);
	}
}

DUMP *unwrap (DUMP *atomEntries, int nAtoms, DUMP simLow, DUMP simHigh)
{
	DUMP *unwrappedCoordinates;
	unwrappedCoordinates = (DUMP *) malloc (nAtoms * sizeof (DUMP));
	float xDim = simHigh.x - simLow.x, yDim = simHigh.y - simLow.y, zDim = simHigh.z - simLow.z;

	for (int i = 0; i < nAtoms; ++i)
	{
		unwrappedCoordinates[i].sino = atomEntries[i].sino;
		unwrappedCoordinates[i].type = atomEntries[i].type;
		unwrappedCoordinates[i].x  = atomEntries[i].x + xDim * (float) (atomEntries[i].ix);
		unwrappedCoordinates[i].y = atomEntries[i].y + yDim * (float) (atomEntries[i].iy);
		unwrappedCoordinates[i].z = atomEntries[i].z + zDim * (float) (atomEntries[i].iz);
		unwrappedCoordinates[i].ix = 0;
		unwrappedCoordinates[i].iy = 0;
		unwrappedCoordinates[i].iz = 0;
	}

	return unwrappedCoordinates;
}

int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		printf("Insufficient args passed. \n\nargv[0] = ./program\nargv[1] = input dump filename\nargv[2] = output dump filename\n\n");
		exit (1);
	}

	FILE *input, *output;
	input = fopen (argv[1], "r");
	output = fopen (argv[2], "w");

	char lineString[1000];
	int lineNumber = 0, nAtoms = 0, currentTimeframe = 0;

	DUMP *atomEntries, com, *recenteredEntries, simLow, simHigh, *unwrappedCoordinates;
	simLow.x = 0; simLow.y = 0; simLow.z = 0;
	simHigh.x = 0; simHigh.y = 0; simHigh.z = 0;
	printf("\n");

	while (fgets (lineString, 1000, input) != NULL)
	{
		lineNumber++;
		currentTimeframe++;

		if (lineNumber <= 8)
			fprintf(output, "%s", lineString);
		if (lineNumber == 9)
			fprintf(output, "ITEM: ATOMS id type x y z ix iy iz\n");

		if (lineNumber == 4 && nAtoms == 0)
		{
			sscanf (lineString, "%d \n", &nAtoms);
			atomEntries = (DUMP *) malloc (nAtoms * sizeof (DUMP));
			recenteredEntries = (DUMP *) malloc (nAtoms * sizeof (DUMP));
			unwrappedCoordinates = (DUMP *) malloc (nAtoms * sizeof (DUMP));
		}

		if (lineNumber == 6 && simLow.x == 0 && simHigh.x == 0)
			sscanf (lineString, "%f %f\n", &simLow.x, &simHigh.x);
		if (lineNumber == 7 && simLow.y == 0 && simHigh.y == 0)
			sscanf (lineString, "%f %f\n", &simLow.y, &simHigh.y);
		if (lineNumber == 8 && simLow.z == 0 && simHigh.z == 0)
			sscanf (lineString, "%f %f\n", &simLow.z, &simHigh.z);

		if (lineNumber > 9 && lineNumber <= nAtoms + 9)
		{
			sscanf (lineString, "%d %d %f %f %f %f %f %f %d %d %d\n", 
				&atomEntries[lineNumber - 10].sino,
				&atomEntries[lineNumber - 10].type,
				&atomEntries[lineNumber - 10].x,
				&atomEntries[lineNumber - 10].y,
				&atomEntries[lineNumber - 10].z,
				&atomEntries[lineNumber - 10].xs,
				&atomEntries[lineNumber - 10].ys,
				&atomEntries[lineNumber - 10].zs,
				&atomEntries[lineNumber - 10].ix,
				&atomEntries[lineNumber - 10].iy,
				&atomEntries[lineNumber - 10].iz);
		}

		if (lineNumber == nAtoms + 9)
		{
			unwrappedCoordinates = unwrap (atomEntries, nAtoms, simLow, simHigh);
			com = findCOM (unwrappedCoordinates, nAtoms);
			recenteredEntries = recenterCoordinates (unwrappedCoordinates, nAtoms, com);
			printCoordinates (recenteredEntries, nAtoms, output);
			lineNumber = 0;

			if ((currentTimeframe % 100) == 0)
			{
				fprintf(stdout, "scanning timeframe: %d...\r", currentTimeframe);
				fflush (stdout);
			}
		}
	}

	fclose (input);
	fclose (output);
	return 0;
}