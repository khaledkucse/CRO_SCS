package DRandIBS_SCS;

import java.io.IOException;
import java.io.RandomAccessFile;


public class DR
{
	int sigSize;
	int n;//number of input strings
	int m;//maxmum lenght of input strings 
	int[] M; //array to keep individual input strings' lenghts
	int[] M_prime;
	class Seq //Sequence
	{
		int[] s;		
		Seq(int pM)
		{
			s= new int[pM];
		}
	}//Seq
	Seq[] S; 	
	String Alphabet = "";	
	int[][] suffix_positions_matrix;
	int retLen; 
	long retTime;
	partSol Global_finalSol ;

	class partSol// Solution
	{
		int[] X;		
		int len;	
		partSol(int pN, int pMaxLen)
		{
			X= new int[pMaxLen];
		}
	}//partSol


	//constructor
	DR    (int Pn,  String pDataFile,  String string_type) throws IOException

	{
		n=Pn;
		int i,j,tmpInx;
		String s;
		m=0; 
		M=new int[n]; 
		S= new Seq[n];
		double MyRand;
		char newChar= ' ';
		if(string_type.equals("D") || string_type.equals("d"))//DNA
			sigSize=4;
		else 
			if(string_type.equals("P") || string_type.equals("p"))//protein
				sigSize=20;

		System.out.println();
		System.out.println("Reading Data from the file: " + pDataFile);
		RandomAccessFile f=new RandomAccessFile(pDataFile, "r");
		f.seek(0);
		i=-1;  
		s = f.readLine();

		while( s != null)
		{
			i++;	
			s.trim();
			M[i]= s.length();
			if(M[i]>m)
				m=M[i];

			S[i]= new Seq(M[i]);

			if(string_type.equals("D") || string_type.equals("d"))//DNA
			{
				for(j=0;j<M[i];j++)
				{
					if(s.charAt(j)=='n' || s.charAt(j)=='N' || s.charAt(j)=='Y' || s.charAt(j)=='R' || s.charAt(j)=='W' || s.charAt(j)=='S' || s.charAt(j)=='K' || s.charAt(j)=='B' || s.charAt(j)=='M')
					{
						if(s.charAt(j)=='Y'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{
								newChar = 'T';
							}else{
								newChar = 'C';
							}
						}
						if(s.charAt(j)=='R'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{
								newChar = 'G';
							}else{
								newChar = 'A';
							}
						}
						if(s.charAt(j)=='W'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{
								newChar = 'A';
							}else{
								newChar = 'T';
							}
						}
						if(s.charAt(j)=='S'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{

								newChar = 'C';
							}else{

								newChar = 'G';
							}
						}
						if(s.charAt(j)=='K'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{
								newChar = 'G';
							}else{
								newChar = 'T';
							}
						}
						if(s.charAt(j)=='M'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{
								newChar = 'A';
							}else{
								newChar = 'C';
							}
						}

						if(s.charAt(j)=='B'){
							MyRand = Math.random();
							if(MyRand < 0.33)
							{
								newChar = 'G';
							}else{
								if(MyRand < 0.66){
									newChar = 'C';
								}
								else
									newChar = 'T';
							}
						}  	

						if(s.charAt(j)=='n' || s.charAt(j)=='N'){
							tmpInx= (int)Math.random()*4;
						}else
							tmpInx= Alphabet.indexOf(newChar);

						if(tmpInx==-1)
						{
							Alphabet = Alphabet + newChar;
							S[i].s[j]= Alphabet.length()-1;
						}
						else
						{
							S[i].s[j]= tmpInx;
						}

					}
					else
					{
						tmpInx= Alphabet.indexOf(s.charAt(j));
						if(tmpInx==-1)
						{
							Alphabet = Alphabet + s.charAt(j);
							S[i].s[j]= Alphabet.length()-1;
						}
						else
						{
							S[i].s[j]= tmpInx;
						}
					}

					if(S[i].s[j]>=sigSize)
					{
						System.out.println("ERROOOOOOOR: unexpected chracter (" + s.charAt(j) + ") encountered in  data: i=" + i + ", j=" +j + " in " + pDataFile);
						System.out.println("Should not happen");
					} 
				} // for j
			} // if DNA seq
			else if(string_type.equals("P") || string_type.equals("p"))//protein
			{
				for(j=0;j<M[i];j++)
				{
					if( s.charAt(j)=='B' || s.charAt(j)=='Z'  || s.charAt(j)=='X'  )

					{
						if(s.charAt(j)=='B'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{

								newChar = 'N';
							}else{

								newChar = 'D';
							}
						}
						if(s.charAt(j)=='Z'){
							MyRand = Math.random();
							if(MyRand < 0.5)
							{

								newChar = 'E';
							}else{

								newChar = 'Q';
							}
						}

						if(s.charAt(j)=='X'){

							tmpInx= (int)Math.random()*20;
						}else
							tmpInx= Alphabet.indexOf(newChar);

						if(tmpInx==-1)
						{
							Alphabet = Alphabet + newChar;
							S[i].s[j]= Alphabet.length()-1;
						}
						else
						{
							S[i].s[j]= tmpInx;
						}
					}
					else
					{
						tmpInx= Alphabet.indexOf(s.charAt(j));
						if(tmpInx==-1)
						{
							Alphabet = Alphabet + s.charAt(j);
							S[i].s[j]= Alphabet.length()-1;
						}
						else
						{
							S[i].s[j]= tmpInx;
						}
					}

					if(S[i].s[j]>=sigSize)
					{
						System.out.println("ERROOOOOOOR: unexpected chracter (" + s.charAt(j) + ") encountered in  data: i=" + i + ", j=" +j + " in " + pDataFile);
						System.out.println("Should not happen");
					}
				} //  for j
			}
			else if (string_type.equals("O")||string_type.equals("o"))//others
			{
				for(j=0;j<M[i];j++)
				{  
					tmpInx= Alphabet.indexOf(s.charAt(j));
					if(tmpInx==-1)
					{
						Alphabet = Alphabet + s.charAt(j);
						S[i].s[j]= Alphabet.length()-1;
					}
					else
					{
						S[i].s[j]= tmpInx;
					}

				}//  for j
				sigSize=Alphabet.length();
			}
			s = f.readLine();
		}


		if(Alphabet.length()!=sigSize)
		{
			System.out.println("ERROOOOOOOR: real sigSize=" + Alphabet.length() + " as oppose to " + sigSize + " in " + pDataFile);			
			return;
		}
		n=i+1;
		f.close();  
	}//constructor
	//-------------------------------------------------

	//-------------------------------------------------
	private void BShPr(partSol pFinalSol)
	{
		retLen=0;	
		int i,i2,i3,j,k,kx,k2,t,counter;

		int[] B_x=new int[sigSize*m];
		int[][] c1_x=new int[sigSize][sigSize*m];
		int[][] c2_x=new int[sigSize*sigSize][sigSize*m];
		int[][] c3_x=new int[sigSize*sigSize*sigSize][sigSize*m];

		int[] B_p=new int[n];
		int[][] c1_p=new int[sigSize][n];
		int[][] c2_p=new int[sigSize*sigSize][n];
		int[][] c3_p=new int[sigSize*sigSize*sigSize][n];

		int[] c3_h=new int[sigSize*sigSize*sigSize];

		int[][] c2_father_p=new int[sigSize*sigSize][n];
		int[][] c3_father_p=new int[sigSize*sigSize*sigSize][n];	

		int nC1,nC2,nC3;//current size of the lists
		int[] max_index=new int[sigSize*sigSize*sigSize];//for the children
		int curLen;
		boolean tmpChildUsable,flgAnswer;	
		int tmpMaxHInx,tmpCount;

		int max_m=M[0];
		for(i=1;i<n;i++)
			if(M[i]>max_m)
				max_m=M[i];


		int[] x_alphabet=new int[sigSize*max_m];
		for(i=0;i<max_m;i++)
			for(i2=0;i<sigSize;i++)
				x_alphabet[i*max_m+i2]=i2;


		for(i2=0;i2<n;i2++)
			B_p[i2]=-1;

		nC3=0;
		curLen=0;
		while(true) //until no more expansion possible
		{	
			counter=0;

			for(t=0;t<sigSize;t++)
			{	
				flgAnswer=true;
				tmpChildUsable=false;
				for(i2=0;i2<n;i2++)
				{

					if( B_p[i2]<M[i2]-1 && S[i2].s[B_p[i2]+1]==t )
					{
						c1_p[counter][i2]= B_p[i2]+1;	
						tmpChildUsable=true;
					}
					else
					{
						c1_p[counter][i2]= B_p[i2];
					}
					if(c1_p[counter][i2] < M[i2]-1)
						flgAnswer=false;
				}
				if(tmpChildUsable)// NOTE the order
				{
					for(kx=0;kx<curLen;kx++)
						c1_x[counter][kx]= B_x[kx];
					c1_x[counter][curLen]= t;//NOTE: starting from 0
					counter++;
				}				
				if(flgAnswer)//prepare the final answer
				{					
					retLen=curLen+1;

					if(retLen>(sigSize*max_m))//alphabet answer is better
					{
						for(j=0;j<=sigSize*max_m;j++)
							pFinalSol.X[j]=x_alphabet[j];
						pFinalSol.len=sigSize*max_m;
					}
					else

					{
						for(j=0;j<=curLen;j++)//Note =
							pFinalSol.X[j]=c1_x[counter-1][j];

						pFinalSol.len=retLen;
					}
					return;
				}
			}

			nC1=counter;

			counter=0;		
			for(k=0;k<nC1;k++)//for each candidate partial solution
			{
				for(t=0;t<sigSize;t++)
				{		
					flgAnswer=true;
					tmpChildUsable=false;
					for(i2=0;i2<n;i2++)
					{
						c2_father_p[counter][i2]=c1_p[k][i2];
						if( c1_p[k][i2]<M[i2]-1 && S[i2].s[c1_p[k][i2]+1]==t )
						{
							c2_p[counter][i2]= c1_p[k][i2]+1;
							tmpChildUsable=true;
						}
						else
						{
							c2_p[counter][i2]= c1_p[k][i2];
						}
						if(c2_p[counter][i2] < M[i2]-1)
							flgAnswer=false;
					}
					if(tmpChildUsable)// NOTE the order
					{
						for(kx=0;kx<curLen+1;kx++)
							c2_x[counter][kx]= c1_x[k][kx];
						c2_x[counter][curLen+1]= t;//NOTE: starting from 0
						counter++;
					}				

					if(flgAnswer)//prepare the final answer
					{					
						retLen=curLen+2;

						if(retLen>(sigSize*max_m))//alphabet answer is better
						{
							for(j=0;j<=sigSize*max_m;j++)
								pFinalSol.X[j]=x_alphabet[j]; 
							pFinalSol.len=sigSize*max_m;
						}
						else

						{
							for(j=0;j<=curLen+1;j++)//Note =
								pFinalSol.X[j]=c2_x[counter-1][j];
							pFinalSol.len=retLen;
						}
						return;
					}				
				}
			}
			nC2=counter;

			counter=0;		
			for(k=0;k<nC2;k++)//for each candidate partial solution
			{
				for(t=0;t<sigSize;t++)
				{		
					flgAnswer=true;
					tmpChildUsable=false;
					for(i2=0;i2<n;i2++)
					{
						c3_father_p[counter][i2]=c2_father_p[k][i2];
						if( c2_p[k][i2]<M[i2]-1 && S[i2].s[c2_p[k][i2]+1]==t )
						{
							c3_p[counter][i2]= c2_p[k][i2]+1;
							tmpChildUsable=true;
						}
						else
						{
							c3_p[counter][i2]= c2_p[k][i2];
						}
						if(c3_p[counter][i2] < M[i2]-1)
							flgAnswer=false;
					}
					if(tmpChildUsable)// NOTE the order
					{
						for(kx=0;kx<curLen+2;kx++)
							c3_x[counter][kx]= c2_x[k][kx];
						c3_x[counter][curLen+2]= t;//NOTE: starting from 0
						counter++;
					}
					if(flgAnswer)//prepare the final answer
					{					
						retLen=curLen+3;

						if(retLen>(sigSize*max_m))//alphabet answer is better
						{
							for(j=0;j<=sigSize*max_m;j++)
								pFinalSol.X[j]=x_alphabet[j]; 
							pFinalSol.len=sigSize*max_m;
						}
						else

						{
							for(j=0;j<=curLen+2;j++)
								pFinalSol.X[j]=c3_x[counter-1][j];
							pFinalSol.len=retLen;
						}
						return;
					}				
				}
			}
			nC3=counter;			


			if(nC3>0)
			{				


				for(k2=0;k2<nC3;k2++)
				{	
					c3_h[k2]=0;
					for(i3=0;i3<n;i3++)
						c3_h[k2]= c3_h[k2] + c3_p[k2][i3];

				}//for k2				  


				tmpCount=0;
				int num_max;


				//initialize tmpMaxHInx to the first available (candidate) solution
				tmpMaxHInx= 0;				
				num_max=0;
				max_index[num_max]=tmpMaxHInx;

				//look for a better one, if any
				for(k2=tmpMaxHInx+1;k2<nC3;k2++)
					if( c3_h[k2] > c3_h[tmpMaxHInx])
					{
						tmpMaxHInx=k2;
						num_max=0;
						max_index[num_max]=k2;
					}
					else
					{
						if( c3_h[k2] == c3_h[tmpMaxHInx])
						{
							num_max++;
							max_index[num_max]=k2;
						}
					}

				tmpMaxHInx=max_index[(int)((num_max+1)* Math.random())];

				if(tmpMaxHInx < nC3)
				{
					//select:
					for(j=0;j<=curLen;j++)//Note =
						B_x[j]=c3_x[tmpMaxHInx][j];
					for(i2=0;i2<n;i2++)
						B_p[i2]=c3_father_p[tmpMaxHInx][i2];

					tmpCount++;
				}
				else
					break;


				nC3=0;
				curLen++;
			}//if nC>0
			else
			{
				System.out.println("Error!!!!!!!!!!!!!!: Unexpected exit point");
				return;
			}//else nC>0
		}//while true
	}//BShPr
	//-------------------------------------------------
	public void mainMethod()
	{
		int i,j,k;

		partSol finalSol= new partSol(n,sigSize*m);
		BShPr(finalSol);

		partSol reduced_sol=new partSol(n, sigSize*m);
		reduced_sol=reduction(finalSol);		
		retLen=reduced_sol.len;						
		finalSol=reduced_sol;	

		// add by ft 26_4
		Global_finalSol = finalSol;

		//test the answer:
		boolean flgOK=true;
		for(i=0; i<n && flgOK; i++)
		{
			k=0;
			for(j=0; j<M[i] && flgOK; j++)
			{
				while(k<finalSol.len && S[i].s[j]!=finalSol.X[k])
					k++;
				if(k>=finalSol.len)
					flgOK=false;
				k++;
			}
		}
		if(!flgOK)
		{
			System.out.println("\n           flgOK=" + flgOK);
			retLen=-1;
		}
	}//mainMethod

	//---------------------------------------------------
	private partSol reduction(partSol sol)
	{
		int i,j,num_of_zeros,reduce=0,counter=0;
		M_prime=new int[n];


		suffix_positions_matrix= new int[n][];
		for(i=0;i<n;i++)
			suffix_positions_matrix[i]=new int[M[i]];

		boolean improvement=true;		
		for(i=0;i<n;i++)		
			M_prime[i]=M[i];			

		while(improvement)
		{
			improvement=false;
			for(i=1;i<sol.len;i=i+1)//i is were the solution is divided intto right and left solutions
			{
				if(i==1)
				{					
					calculate_suffix_matrix_2(sol,M,reduce);						

					M=longest_suffix_first_time();
				}
				else
					M=longest_suffix_2(sol, i,M);		


				num_of_zeros=0;
				for(j=0;j<n;j++)
					if(M[j]==0)
						num_of_zeros++;

				partSol left_sol=new partSol(n, sigSize*m);
				if(num_of_zeros==n)
				{
					left_sol.len=0;
				}
				else
				{	
					BShPr( left_sol);
				}
				if(left_sol.len < i)
				{
					reduce=i-left_sol.len;
					counter++;
					improvement=true;
					for(j=0;j<left_sol.len;j++)
						sol.X[j]=left_sol.X[j];
					for(j=i;j<sol.len;j++)
						sol.X[j+left_sol.len-i]=sol.X[j];
					sol.len=sol.len - i + left_sol.len;	
					break;
				}							
			}
		}
		for(j=0;j<n;j++)		
			M[j]=M_prime[j];

		return sol;
	}

	private int[] longest_suffix_2(partSol sol,int half , int[] previous_suffix)
	{
		int i;    	   	
		int[] l=new int[n];
		for(i=0;i<n;i++)
		{
			if(previous_suffix[i]==M_prime[i])
				l[i]=previous_suffix[i];
			else
			{
				if(half -1 ==suffix_positions_matrix[i][previous_suffix[i]])	
					l[i]=previous_suffix[i]+1; //suffix of String[i] is :  s[l[i]]s[l[i]+1]...s[len-1]
				else
					l[i]=previous_suffix[i];   //suffix of String[i] is :  s[l[i]]s[l[i]+1]...s[len-1] 
			}
		}
		return l;
	}

	private int[] longest_suffix_first_time()
	{
		int i;    	   	
		int[] l=new int[n];
		for(i=0;i<n;i++)
		{
			if(suffix_positions_matrix[i][0]==0)
				l[i]=1;
			else
				l[i]=0;   		
		}
		return l;
	}

	private void calculate_suffix_matrix_2(partSol sol,int[] previous_suffix,int reduction)
	{
		int i,j,pos;    	
		for(i=0;i<n;i++)
		{    	
			pos=sol.len -1;  
			for(j=M_prime[i]-1;j>=previous_suffix[i];j--)
			{
				suffix_positions_matrix[i][j]=suffix_positions_matrix[i][j]-reduction;
				pos = suffix_positions_matrix[i][j] -1;
			}    			

			for(j=previous_suffix[i]-1;j>=0;j--)
			{
				while (S[i].s[j]!=sol.X[pos]) 
				{
					pos--;
					if(pos==-1)
						break;
				} 
				if(pos==-1)
				{
					break;
				}	
				suffix_positions_matrix[i][j]=pos;//character num j of string num i is supported by position "pos" of solution, if searching from right
				pos--;	
				if(pos==-1)
				{
					break;
				}	
			} 
		} 
	}
}//DR
