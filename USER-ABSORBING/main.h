#include <cmath>
#include <iostream>

using namespace std;
int width = 5;
int heigh = 5;
int deep = 5;
double lattice = 1.;
double box[3] = {width*lattice, heigh*sqrt(3.)*lattice};

void Create_Side(int width, int , vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize(width);
	innerPos.resize(width);
	for(int i=0; i<width; ++i)
	{
		bcPos[i][0] = i;
		bcPos[i][1] = 0.;
		bcPos[i][2] = 0.;

		innerPos[i][0] = i+0.5;
		innerPos[i][1] = sqrt(3.)*0.5;
		innerPos[i][2] = 0.;
	}
}


void Create_SideX_FCC(int, int heigh, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize( heigh*deep*4 );
	innerPos.resize( heigh*deep*4 );
	for(int i=0; i<heigh; ++i)
	for(int j=0; j<deep; ++j)
	{
		bcPos[i*deep*4+j*4] = {0., i, j};
		bcPos[i*deep*4+j*4+1] = {0.5, i+0.5, j};
		bcPos[i*deep*4+j*4+2] = {0.5, i, j+0.5};
		bcPos[i*deep*4+j*4+3] = {0., i+0.5, j+0.5};

		innerPos[i*deep*4+j*4] = {1., i, j};
		innerPos[i*deep*4+j*4+1] = {1.5, i+0.5, j};
		innerPos[i*deep*4+j*4+2] = {1.5, i, j+0.5};
		innerPos[i*deep*4+j*4+3] = {1., i+0.5, j+0.5};
	}
}

void Create_SideY_Si(int width, int, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize( width*deep*8 );
	innerPos.resize( width*deep*8 );
	for ( int i=0; i<width; ++i )
	for ( int j=0; j<deep; ++j )
	{
		bcPos[i*deep*8+j*8] = {i, 0., j};
		bcPos[i*deep*8+j*8+1] = {i+0.5, 0.5, j};
		bcPos[i*deep*8+j*8+2] = {i+0.5, 0., j+0.5};
		bcPos[i*deep*8+j*8+3] = {i, 0.5, j+0.5};
		bcPos[i*deep*8+j*8+4] = {i+0.25, 0.25, j+0.25};
		bcPos[i*deep*8+j*8+5] = {i+0.75, 0.75, j+0.25};
		bcPos[i*deep*8+j*8+6] = {i+0.25, 0.75, j+0.75};
		bcPos[i*deep*8+j*8+7] = {i+0.75, 0.25, j+0.75};

		innerPos[i*deep*8+j*8] = {i, -1., j};
		innerPos[i*deep*8+j*8+1] = {i+0.5, -0.5, j};
		innerPos[i*deep*8+j*8+2] = {i+0.5, -1., j+0.5};
		innerPos[i*deep*8+j*8+3] = {i, -0.5, j+0.5};
		innerPos[i*deep*8+j*8+4] = {i+0.25, -0.75, j+0.25};
		innerPos[i*deep*8+j*8+5] = {i+0.75, -0.25, j+0.25};
		innerPos[i*deep*8+j*8+6] = {i+0.25, -0.25, j+0.75};
		innerPos[i*deep*8+j*8+7] = {i+0.75, -0.75, j+0.75};
/*		innerPos[i*deep*8+j*8] = {i, 1., j};
		innerPos[i*deep*8+j*8+1] = {i+0.5, 1.5, j};
		innerPos[i*deep*8+j*8+2] = {i+0.5, 1., j+0.5};
		innerPos[i*deep*8+j*8+3] = {i, 1.5, j+0.5};
		innerPos[i*deep*8+j*8+4] = {i+0.25, 1.25, j+0.25};
		innerPos[i*deep*8+j*8+5] = {i+0.75, 1.75, j+0.25};
		innerPos[i*deep*8+j*8+6] = {i+0.25, 1.75, j+0.75};
		innerPos[i*deep*8+j*8+7] = {i+0.75, 1.25, j+0.75};*/
	}
}


void Create_SideZ_Si(int width, int height, int, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
	bcPos.resize( width*height*8 );
	innerPos.resize( width*height*8 );
	for (int i=0; i<width; ++i)
	for (int j=0; j<height; ++j)
	{
		bcPos[i*height*8+j*8] = {i, j, 0.};
		bcPos[i*height*8+j*8+1] = {i+0.5, j, 0.5};
		bcPos[i*height*8+j*8+2] = {i+0.5, j+0.5, 0.};
		bcPos[i*height*8+j*8+3] = {i, j+0.5, 0.5};
		bcPos[i*height*8+j*8+4] = {i+0.25, j+0.25, 0.25};
		bcPos[i*height*8+j*8+5] = {i+0.75, j+0.25, 0.75};
		bcPos[i*height*8+j*8+6] = {i+0.25, j+0.75, 0.75};
		bcPos[i*height*8+j*8+7] = {i+0.75, j+0.75, 0.25};

		innerPos[i*height*8+j*8] = {i, j, -1.};
		innerPos[i*height*8+j*8+1] = {i+0.5, j, -0.5};
		innerPos[i*height*8+j*8+2] = {i+0.5, j+0.5, -1.};
		innerPos[i*height*8+j*8+3] = {i, j+0.5, -0.5};
		innerPos[i*height*8+j*8+4] = {i+0.25, j+0.25, -0.75};
		innerPos[i*height*8+j*8+5] = {i+0.75, j+0.25, -0.25};
		innerPos[i*height*8+j*8+6] = {i+0.25, j+0.75, -0.25};
		innerPos[i*height*8+j*8+7] = {i+0.75, j+0.75, -0.75};
/*		innerPos[i*height*8+j*8] = {i, j, 1.};
		innerPos[i*height*8+j*8+1] = {i+0.5, j, 1.5};
		innerPos[i*height*8+j*8+2] = {i+0.5, j+0.5, 1.};
		innerPos[i*height*8+j*8+3] = {i, j+0.5, 1.5};
		innerPos[i*height*8+j*8+4] = {i+0.25, j+0.25, 1.25};
		innerPos[i*height*8+j*8+5] = {i+0.75, j+0.25, 1.75};
		innerPos[i*height*8+j*8+6] = {i+0.25, j+0.75, 1.75};
		innerPos[i*height*8+j*8+7] = {i+0.75, j+0.75, 1.25};*/
	}
}


void Create_SideX_Si(int, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize( height*deep*8 );
	innerPos.resize( height*deep*8 );
	for (int i=0; i<height; ++i)
	for (int j=0; j<deep; ++j)
	{
		bcPos[i*deep*8+j*8] = {0., i, j};
		bcPos[i*deep*8+j*8+1] = {0.5, i+0.5, j};
		bcPos[i*deep*8+j*8+2] = {0., i+0.5, j+0.5};
		bcPos[i*deep*8+j*8+3] = {0.5, i, j+0.5};
		bcPos[i*deep*8+j*8+4] = {0.25, i+0.25, j+0.25};
		bcPos[i*deep*8+j*8+5] = {0.75, i+0.75, j+0.25};
		bcPos[i*deep*8+j*8+6] = {0.75, i+0.25, j+0.75};
		bcPos[i*deep*8+j*8+7] = {0.25, i+0.75, j+0.75};

		innerPos[i*deep*8+j*8] = {-1., i, j};
		innerPos[i*deep*8+j*8+1] = {-0.5, i+0.5, j};
		innerPos[i*deep*8+j*8+2] = {-1., i+0.5, j+0.5};
		innerPos[i*deep*8+j*8+3] = {-0.5, i, j+0.5};
		innerPos[i*deep*8+j*8+4] = {-0.75, i+0.25, j+0.25};
		innerPos[i*deep*8+j*8+5] = {-0.25, i+0.75, j+0.25};
		innerPos[i*deep*8+j*8+6] = {-0.25, i+0.25, j+0.75};
		innerPos[i*deep*8+j*8+7] = {-0.75, i+0.75, j+0.75};
/*		innerPos[i*deep*8+j*8] = {1., i, j};
		innerPos[i*deep*8+j*8+1] = {1.5, i+0.5, j};
		innerPos[i*deep*8+j*8+2] = {1., i+0.5, j+0.5};
		innerPos[i*deep*8+j*8+3] = {1.5, i, j+0.5};
		innerPos[i*deep*8+j*8+4] = {1.25, i+0.25, j+0.25};
		innerPos[i*deep*8+j*8+5] = {1.75, i+0.75, j+0.25};
		innerPos[i*deep*8+j*8+6] = {1.75, i+0.25, j+0.75};
		innerPos[i*deep*8+j*8+7] = {1.25, i+0.75, j+0.75};*/
	}
}


void Create_EdgeXY_Si ( int width, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
	bcPos.resize( deep*(2*width+2*height-4)*8 );
	innerPos.resize( deep*(2*width+2*height-12)*8 );
	int index = 0;
	for ( int i=0; i<deep; ++i )
	for ( int j=0; j<height; ++j )
	for ( int k=0; k<width; ++k )
	{
		if (j>0 && j<height-1 && k>0 && k<width-1) continue;
		bcPos[index++] = {k, j, i};
		bcPos[index++] = {0.5+k, 0.5+j, i};
		bcPos[index++] = {0.5+k, j, 0.5+i};
		bcPos[index++] = {k, 0.5+j, 0.5+i};
		bcPos[index++] = {0.25+k, 0.25+j, 0.25+i};
		bcPos[index++] = {0.75+k, 0.75+j, 0.25+i};
		bcPos[index++] = {0.25+k, 0.75+j, 0.75+i};
		bcPos[index++] = {0.75+k, 0.25+j, 0.75+i};
	}

	index = 0;
	for ( int i=0; i<deep; ++i )
	for ( int j=1; j<height-1; ++j )
	for ( int k=1; k<width-1; ++k )
	{
		if (j>1 && j<height-2 && k>1 && k<width-2) continue;
		innerPos[index++] = {k, j, i};
		innerPos[index++] = {0.5+k, 0.5+j, i};
		innerPos[index++] = {0.5+k, j, 0.5+i};
		innerPos[index++] = {k, 0.5+j, 0.5+i};
		innerPos[index++] = {0.25+k, 0.25+j, 0.25+i};
		innerPos[index++] = {0.75+k, 0.75+j, 0.25+i};
		innerPos[index++] = {0.25+k, 0.75+j, 0.75+i};
		innerPos[index++] = {0.75+k, 0.25+j, 0.75+i};
	}
}


void Create_EdgeXZ_Si ( int width, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
	bcPos.resize( height*(2*width+2*deep-4)*8 );
	innerPos.resize( height*(2*width+2*deep-12)*8 );
	int index = 0;
	for ( int i=0; i<height; ++i )
	for ( int j=0; j<deep; ++j )
	for ( int k=0; k<width; ++k )
	{
		if (j>0 && j<height-1 && k>0 && k<width-1) continue;
		bcPos[index++] = {k, i, j};
		bcPos[index++] = {0.5+k, 0.5+i, j};
		bcPos[index++] = {0.5+k, i, 0.5+j};
		bcPos[index++] = {k, 0.5+i, 0.5+j};
		bcPos[index++] = {0.25+k, 0.25+i, 0.25+j};
		bcPos[index++] = {0.75+k, 0.75+i, 0.25+j};
		bcPos[index++] = {0.25+k, 0.75+i, 0.75+j};
		bcPos[index++] = {0.75+k, 0.25+i, 0.75+j};
	}

	index = 0;
	for ( int i=0; i<height; ++i )
	for ( int j=1; j<deep-1; ++j )
	for ( int k=1; k<width-1; ++k )
	{
		if (j>1 && j<height-2 && k>1 && k<width-2) continue;
		innerPos[index++] = {k, i, j};
		innerPos[index++] = {0.5+k, 0.5+i, j};
		innerPos[index++] = {0.5+k, i, 0.5+j};
		innerPos[index++] = {k, 0.5+i, 0.5+j};
		innerPos[index++] = {0.25+k, 0.25+i, 0.25+j};
		innerPos[index++] = {0.75+k, 0.75+i, 0.25+j};
		innerPos[index++] = {0.25+k, 0.75+i, 0.75+j};
		innerPos[index++] = {0.75+k, 0.25+i, 0.75+j};
	}
}


void Create_EdgeYZ_Si ( int width, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
	bcPos.resize( width*(2*height+2*deep-4)*8 );
	innerPos.resize( width*(2*height+2*deep-12)*8 );
	int index = 0;
	for ( int i=0; i<width; ++i )
	for ( int j=0; j<deep; ++j )
	for ( int k=0; k<height; ++k )
	{
		if (j>0 && j<height-1 && k>0 && k<width-1) continue;
		bcPos[index++] = {i, k, j};
		bcPos[index++] = {0.5+i, 0.5+k, j};
		bcPos[index++] = {0.5+i, k, 0.5+j};
		bcPos[index++] = {i, 0.5+k, 0.5+j};
		bcPos[index++] = {0.25+i, 0.25+k, 0.25+j};
		bcPos[index++] = {0.75+i, 0.75+k, 0.25+j};
		bcPos[index++] = {0.25+i, 0.75+k, 0.75+j};
		bcPos[index++] = {0.75+i, 0.25+k, 0.75+j};
	}

	index = 0;
	for ( int i=0; i<width; ++i )
	for ( int j=1; j<deep-1; ++j )
	for ( int k=1; k<height-1; ++k )
	{
		if (j>1 && j<height-2 && k>1 && k<width-2) continue;
		innerPos[index++] = {i, k, j};
		innerPos[index++] = {0.5+i, 0.5+k, j};
		innerPos[index++] = {0.5+i, k, 0.5+j};
		innerPos[index++] = {i, 0.5+k, 0.5+j};
		innerPos[index++] = {0.25+i, 0.25+k, 0.25+j};
		innerPos[index++] = {0.75+i, 0.75+k, 0.25+j};
		innerPos[index++] = {0.25+i, 0.75+k, 0.75+j};
		innerPos[index++] = {0.75+i, 0.25+k, 0.75+j};
	}
}


void Create_Corner_Si ( int width, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
	bcPos.resize( (width*height*2 + (deep-2)*(width*2+height*2-4))*8 );
	innerPos.resize( ((width-2)*(height-2)*2 + (deep-4)*((width-2)*2+(height-2)*2-4))*8 );

	int index = 0;
	for ( int i=0; i<width; ++i )
	for ( int j=0; j<height; ++j )
	for ( int k=0; k<deep; ++k )
	{
		if ( i>0 && i<width-1 && j>0 && j<height-1 && k>0 && k<deep-1 ) continue;
		bcPos[index++] = {i, j, k};
		bcPos[index++] = {0.5+i, 0.5+j, k};
		bcPos[index++] = {0.5+i, j, 0.5+k};
		bcPos[index++] = {i, 0.5+j, 0.5+k};
		bcPos[index++] = {0.25+i, 0.25+j, 0.25+k};
		bcPos[index++] = {0.75+i, 0.75+j, 0.25+k};
		bcPos[index++] = {0.25+i, 0.75+j, 0.75+k};
		bcPos[index++] = {0.75+i, 0.25+j, 0.75+k};
	}

	index = 0;
	for ( int i=1; i<width-1; ++i )
	for ( int j=1; j<height-1; ++j )
	for ( int k=1; k<deep-1; ++k )
	{
		if ( i>1 && i<width-2 && j>1 && j<height-2 && k>1 && k<deep-2 ) continue;
		innerPos[index++] = {i, j, k};
		innerPos[index++] = {0.5+i, 0.5+j, k};
		innerPos[index++] = {0.5+i, j, 0.5+k};
		innerPos[index++] = {i, 0.5+j, 0.5+k};
		innerPos[index++] = {0.25+i, 0.25+j, 0.25+k};
		innerPos[index++] = {0.75+i, 0.75+j, 0.25+k};
		innerPos[index++] = {0.25+i, 0.75+j, 0.75+k};
		innerPos[index++] = {0.75+i, 0.25+j, 0.75+k};
	}
}


void Create_CornerP_Si ( int width, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos )
{
//	bcPos.resize( (width*height*2 + (deep-2)*(width*2+height*2-4))*8 );
//	innerPos.resize( ((width-2)*(height-2)*2 + (deep-4)*((width-2)*2+(height-2)*2-4))*8 );
//	int index = 0;
}


void Create_SideV(int , int heigh, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize(heigh*2);
	innerPos.resize(heigh*2);
	for(int i=0; i<heigh; ++i)
	{
		bcPos[2*i][0] = 0.;
		bcPos[2*i][1] = i*sqrt(3.);
		bcPos[2*i][2] = 0.;
		bcPos[2*i+1][0] = 0.5;
		bcPos[2*i+1][1] = (i+0.5)*sqrt(3.);
		bcPos[2*i+1][2] = 0.;

		innerPos[2*i][0] = 1.;
		innerPos[2*i][1] = i*sqrt(3.);
		innerPos[2*i][2] = 0.;
		innerPos[2*i+1][0] = 1.5;
		innerPos[2*i+1][1] = (i+0.5)*sqrt(3.);
		innerPos[2*i+1][2] = 0.;
	}
}


void Create_SideV_Si(int, int height, int deep, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	bcPos.resize(height*8);
	innerPos.resize(height*8);
	for (int i=0; i<height; ++i)
	for (int j=0; j<deep; ++j)
	{
		bcPos[i*deep*8+j*8] = {0., i, j};
		bcPos[i*deep*8+j*8+1] = {0.5, i+0.5, j};
		bcPos[i*deep*8+j*8+2] = {0.5, i, 0.5+j};
		bcPos[i*deep*8+j*8+3] = {0., i+0.5, 0.5+j};
		bcPos[i*deep*8+j*8+4] = {0.25, i+0.25, 0.25+j};
		bcPos[i*deep*8+j*8+5] = {0.75, i+0.75, 0.25+j};
		bcPos[i*deep*8+j*8+6] = {0.25, i+0.75, 0.75+j};
		bcPos[i*deep*8+j*8+7] = {0.75, i+0.25, 0.75+j};

		innerPos[i*deep*8+j*8] = {1., i, j};
		innerPos[i*deep*8+j*8+1] = {1.5, i+0.5, j};
		innerPos[i*deep*8+j*8+2] = {1.5, i, 0.5+j};
		innerPos[i*deep*8+j*8+3] = {1., i+0.5, 0.5+j};
		innerPos[i*deep*8+j*8+4] = {1.25, i+0.25, 0.25+j};
		innerPos[i*deep*8+j*8+5] = {1.75, i+0.75, 0.25+j};
		innerPos[i*deep*8+j*8+6] = {1.25, i+0.75, 0.75+j};
		innerPos[i*deep*8+j*8+7] = {1.75, i+0.25, 0.75+j};
	}
}


void Create_CornerLO(int width, int heigh, vector<data3d>& pos)
{
	pos.resize(2*width+4*heigh-4);
	int index = 0;
	for(int i=0; i<heigh; ++i)
	{
		pos[index][0] = 0.;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		pos[index][0] = 0.5;
		pos[index][1] = (0.5+i)*sqrt(3.);
		pos[index][2] = 0.;
		++index;

		if(i==0)
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j;
			pos[index][1] = sqrt(3.)*i;
			pos[index][2] = 0.;
			++index;
		}
		if(i==heigh-1)
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j+0.5;
			pos[index][1] = sqrt(3.)*(i+0.5);
			pos[index][2] = 0.;
			++index;
		}

		pos[index][0] = width-1.;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		pos[index][0] = width-0.5;
		pos[index][1] = (0.5+i)*sqrt(3.);
		pos[index][2] = 0.;
		++index;
	}
}


void Create_CornerLL(int width, int heigh, vector<data3d>& pos)
{
	pos.resize(2*width+4*heigh-6);
	int index = 0;
	for(int i=0; i<heigh; ++i)
	{
		pos[index][0] = 0.;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		if(i != heigh-1)
		{
			pos[index][0] = 0.5;
			pos[index][1] = (i+0.5)*sqrt(3.);
			pos[index][2] = 0.;
			++index;
		}

		//--------
		if(i==0)
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j;
			pos[index][1] = 0.;
			pos[index][2] = 0.;
			++index;
		}

		if(i==heigh-1)
		{
		--index;
		for(int j=0; j<width-1; ++j)
		{
			pos[index][0] = j;
			pos[index][1] = sqrt(3.)*i;
			pos[index][2] = 0.;
			++index;
		}
		}
		//----------

		pos[index][0] = width-1.;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		if(i != heigh-1)
		{ 
			pos[index][0] = width-1.5;
			pos[index][1] = (0.5+i)*sqrt(3.);
			pos[index][2] = 0.;
			++index;
		}
	}
}

void Create_CornerOL(int width, int heigh, vector<data3d>& pos)
{
	pos.resize(2*width+4*heigh-4);
	int index = 0;
	for(int i=0; i<heigh; ++i)
	{
		pos[index][0] = 0.5;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		pos[index][0] = 0.;
		pos[index][1] = (i+0.5)*sqrt(3.);
		pos[index][2] = 0.;
		++index;

		if(i==0)
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j+0.5;
			pos[index][1] = 0.;
			pos[index][2] = 0.;
			++index;
		}
		if(i==heigh-1)
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j;
			pos[index][1] = (heigh-0.5)*sqrt(3.);
			pos[index][2] = 0.;
			++index;
		}

		pos[index][0] = width - 0.5;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		pos[index][0] = width-1.;
		pos[index][1] = (i+0.5)*sqrt(3.);
		pos[index][2] = 0.;
		++index;
	}
}


void Create_CornerOO(int width, int heigh, vector<data3d>& pos)
{
	pos.resize(2*width+4*heigh-8);
	int index = 0;
	for(int i=0; i<heigh; ++i)
	{
		pos[index][0] = 0.5;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		if(i != heigh-1)
		{
			pos[index][0] = 0;
			pos[index][1] = (i+0.5)*sqrt(3.);
			pos[index][2] = 0.;
			++index;
		}

		//--------
		if(i==0)
		for(int j=1; j<width-2; ++j)
		{
			pos[index][0] = j+0.5;
			pos[index][1] = 0.;
			pos[index][2] = 0.;
			++index;
		}

		if(i==heigh-1)
		{
		--index;
		for(int j=1; j<width-1; ++j)
		{
			pos[index][0] = j-0.5;
			pos[index][1] = sqrt(3.)*i;
			pos[index][2] = 0.;
			++index;
		}
		}
		//----------

		pos[index][0] = width-1.5;
		pos[index][1] = i*sqrt(3.);
		pos[index][2] = 0.;
		++index;
		if(i != heigh-1)
		{ 
			pos[index][0] = width-1.;
			pos[index][1] = (0.5+i)*sqrt(3.);
			pos[index][2] = 0.;
			++index;
		}
	}
}


void Create_Corner1(int width, int heigh, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	Create_CornerLO(width, heigh, bcPos);
	Create_CornerOL(width-2, heigh-1, innerPos);
	data3d shift;
	shift[0] = 1.;
	shift[1] = sqrt(3.)*0.5;
	for(unsigned i=0; i<innerPos.size(); ++i)
	{
		innerPos[i][0] += shift[0];
		innerPos[i][1] += shift[1];
	}
}


void Create_Corner1_(int width, int heigh, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	Create_CornerOO(width, heigh, bcPos);
	Create_CornerLL(width-2, heigh-1, innerPos);
	data3d shift;
	shift[0] = 1.;
	shift[1] = sqrt(3.)*0.5;
	for(unsigned i=0; i<innerPos.size(); ++i)
	{
		innerPos[i][0] += shift[0];
		innerPos[i][1] += shift[1];
	}
}


void Create_Corner2(int width, int heigh, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	Create_CornerOL(width, heigh, bcPos);
	Create_CornerLO(width-2, heigh-1, innerPos);
	data3d shift;
	shift[0] = 1.;
	shift[1] = sqrt(3.)*0.5;
	for(unsigned i=0; i<innerPos.size(); ++i)
	{
		innerPos[i][0] += shift[0];
		innerPos[i][1] += shift[1];
	}
}


void Create_Corner2_(int width, int heigh, vector<data3d>& bcPos, vector<data3d>& innerPos)
{
	Create_CornerLL(width, heigh, bcPos);
	Create_CornerOO(width-2, heigh-1, innerPos);
	data3d shift;
	shift[0] = 1.;
	shift[1] = sqrt(3.)*0.5;
	for(unsigned i=0; i<innerPos.size(); ++i)
	{
		innerPos[i][0] += shift[0];
		innerPos[i][1] += shift[1];
	}
}


struct PBC : public PBCDifference {
	void
	Correction ( double* diff ) const noexcept override
	{
		if(diff[0] > box[0]*0.5) diff[0]-=box[0];
		if(diff[0] < -box[0]*0.5) diff[0]+=box[0];
		if(diff[1] > box[1]*0.5) diff[1]-=box[1];
		if(diff[1] < -box[1]*0.5) diff[1]+=box[1];
		if(diff[2] > box[2]*0.5) diff[2]-=box[2];
		if(diff[2] < -box[2]*0.5) diff[2]+=box[2];
	}
};


struct NBC : public PBCDifference {
	void
	Correction ( double* ) const noexcept override
	{}
};

