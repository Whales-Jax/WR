/*
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;
int dummy[50];
int dummy_array(int K, int Ntubes);

int main()
{
	int i, j, K, x, y, Njunc, Nrow, Ncol, Ntubes, Nloop = 0;
	int tube_mtrx[100][100];
	int c_pos[50], c_neg[50], pos[50], neg[50], sumrow[50];
	int dummy_sum = 0;
	bool search, forward, reverse;

	//Tube-tube connectivity matrix (TTM)//
	cout << "Enter the number of rows and columns, separated by space." << endl;
	cin >> Nrow >> Ncol;

	Ntubes = Nrow * Ncol;
	Njunc = Ntubes - 1;						//no splits, no merges

	//TTM initialization
	for (i = 1; i <= Ntubes; i++)
	{
		dummy[i] = 0;

		for (j = 1; j <= Ntubes; j++)
		{
			tube_mtrx[i][j] = 0;
		}
	}

	//Input of Tube Connectivity
	for (i = 1; i <= Njunc; i++)
	{
		cout << "\nEnter the tubes connected (Note: source destination)" << endl;
		cin >> x >> y;

		tube_mtrx[x][y] = 1;
		tube_mtrx[y][x] = -1;
	}

	cout << "\n";

	//Search for inlet and outlet tubes, and splitting and merging
	for (i = 1; i <= Ntubes; i++)
	{
		c_pos[i] = 0;
		c_neg[i] = 0;
		pos[i] = 0;
		neg[i] = 0;

		for (j = 1; j <= Ntubes; j++)
		{
			if (tube_mtrx[i][j] == 1)
			{
				c_pos[i]++;
				pos[i]++;
			}

			if (tube_mtrx[i][j] == -1)
			{
				c_neg[i]++;
				neg[i]++;
			}
		}

		if (c_neg[i] == 0) cout << "The input tube is at tube [" << i << "]." << endl;
		if (c_pos[i] == 0) cout << "The output tube is at tube [" << i << "]." << endl;

		if (c_neg[i] > 1) cout << "Merging at tube [" << i << "]." << endl;
		if (c_pos[i] > 1) cout << "Splitting at tube [" << i << "]." << endl;
	}

	cout << "\n\n";

	K = 0;
	search = true;
	forward = false;
	reverse = false;

	do
	{	
		dummy_sum = Ntubes;

		if (search == true)
		{
			Nloop++;

			for (i = 1; i <= Ntubes; i++)
			{
				bool check;
				check = false;

				for (j = 1; j <= Ntubes; j++)
				{
					if (dummy[j] == i)
					{
						check = true;
						break;
					}

					else {}
				}
				
				if (check == false)
				{
					K = i;
					dummy_array(K, Ntubes);
					break;
				}
			}

			search = false;

			if (pos[K] > 0)
			{
				forward = true;
			}

			else
			{
				reverse = true;
			}
		}

		if (forward == true)
		{
			//cout << "\nHello F1";
			if (pos[K] == 0 || c_pos[K] == 0)			//All of the branches of the split have been tranversed OR K is the outlet tube
			{
				//cout << "\nHello F2";
				if (c_pos[K] == 0 && neg[K] > 0)
				{
					//cout << "\nHello F3";
					forward = false;
					reverse = true;
				}

				else
				{
					//cout << "\nHello F4";
					for (j = Ntubes; j >= 1; j--)			//Looks for a new path from the dummy array
					{
						bool other;
						other = true;

						if (dummy[j] != 0)
						{
							//cout << "\nHello F5";
							x = dummy[j];

							if (c_pos[x] > 1 && sumrow[x] != 0)
							{
								//cout << "\nHello F6";
								other = false;
								K = x;
								break;
							}

							if (c_neg[x] > 1 || sumrow[x] != 0)
							{
								//cout << "\nHello F7";
								other = false;
								K = x;
								forward = false;
								reverse = true;
								break;
							}

						}

						if (j == 1 && other == true)
						{
							//cout << "\nHello F8";
							forward = false;
							search = true;
						}

						else {}
					}
				}
			}

			else
			{
				//cout << "\nHello F9";
				for (i = 1; i <= Ntubes; i++)
				{
					if (tube_mtrx[K][i] == 1)
					{
						//cout << "\nHello F10";
						if (c_pos[K] > 1)
						{
							//cout << "\nHello F11";
							bool check;
							check = false;

							for (j = 1; j <= Ntubes; j++)
							{
								if (dummy[j] == i)
								{
									//cout << "\nHello F12";
									check = true;
									break;
								}

								else {}
							}

							if (check == false)
							{
								//cout << "\nHello F13";
								pos[K]--;
								neg[i]--;
								K = i;
								dummy_array(K, Ntubes);
								break;
							}
							
						}

						else
						{
							//cout << "\nHello F14";
							dummy_array(i, Ntubes);
							pos[K]--;
							neg[i]--;
							K = i;
							break;
						}
					}

					else {}
				}
			}
		}

		if (reverse == true)
		{
			//cout << "\nHello R0";
			if (neg[K] == 0 || c_neg[K] == 0)			//All of the branches of the split have been tranversed OR K is the outlet tube
			{
				//cout << "\nHello R1";
				if (c_neg[K] == 0 && pos[K] > 0)
				{
					//cout << "\nHello R2";
					reverse = false;
					forward = true;
				}
				
				else
				{
					//cout << "\nHello R3";
					for (j = Ntubes; j >= 1; j--)			//Looks for a new path from the dummy array
					{
						bool other;
						other = true;

						if (dummy[j] != 0)
						{
							//cout << "\nHello R4";
							x = dummy[j];

							if (c_neg[x] > 1 && sumrow[x] != 0)						//before  x != K
							{
								//cout << "\nHello R5";
								other = false;
								K = x;
								break;
							}

							if (c_pos[x] > 1 && sumrow[x] != 0)
							{
								//cout << "\nHello R6";
								other = false;
								K = x;
								reverse = false;
								forward = true;
								break;
							}
						}

						if (j == 1 && other == true)
						{
							//cout << "\nHello R7";
							reverse = false;
							search = true;
						}

						else {}
					}
				}
			}

			else
			{
				//cout << "\nHello R8";
				for (i = 1; i <= Ntubes; i++)
				{
					if (tube_mtrx[K][i] == -1)
					{
						//cout << "\nHello R9";
						if (c_neg[K] > 1)
						{
							//cout << "\nHello R10";
							bool check;
							check = false;

							for (j = 1; j <= Ntubes; j++)
							{
								if (dummy[j] == i)
								{
									//cout << "\nHello R11";
									check = true;
									break;
								}

								else {}
							}

							if (check == false)
							{
								//cout << "\nHello R12";
								neg[K]--;
								pos[i]--;
								K = i;
								dummy_array(K, Ntubes);
								break;
							}

						}

						else
						{
							//cout << "\nHello R13";
							dummy_array(i, Ntubes);
							neg[K]--;
							pos[i]--;
							K = i;
							break;
						}
					}

					else {}
				}
			}
		}

		//Dummy_sum counter											//Counts how many unique tubes have been traversed
		for (i = 1; i <= Ntubes; i++)
		{
			if (dummy[i] != 0)
			{
				dummy_sum--;
			}

			else {}
		}


		//cout << "\nDUMMY\n";
		for (i = 1; i <= Ntubes; i++)
		{
			sumrow[i] = pos[i] + neg[i];
			//cout << "\t" << dummy[i];
		}

			   
	} while (dummy_sum > 0);

	cout << "\n\nDUMMY\t";

	for (i = 1; i <= Ntubes; i++)
	{
		cout << "\t" << dummy[i];
	}

	cout << "\n\nNumber of parallel circuit(s) = " << Nloop << "\n";

	cout << "\n\nCount";
	for (i = 1; i <= Ntubes; i++)
	{
		cout << "\nTube " << i << "\tpos = " << pos[i] << "\tneg = " << neg[i] << "\tsum row = " << sumrow[i];
	}

	cout << "\n\n";
	system("pause");
	return 0;
}

int dummy_array(int K, int Ntubes)									//Allocates a value into the dummy array
{
	int i, j;
	bool check;

	check = false;

	for (i = 1; i <= Ntubes; i++)
	{
		if (dummy[i] == 0)											//Puts K in an "empty" slot of the dummy array
		{
			for (j = 1; j <= i; j++)								//Checks if K has already been in the dummy array
			{
				if (dummy[j] == K)
				{
					check = true;
					break;
				}

				else {}
			}

			if (check == false)
			{
				dummy[i] = K;
			}

			else {}

			break;
		}

		else {}
	}
	return 0;
}*/