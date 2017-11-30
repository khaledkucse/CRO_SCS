package DRandIBS_SCS;

import java.io.IOException;
import java.io.RandomAccessFile;



public class IBS_SCS
{
	int sigSize;
	int n;//number of input strings
	int m;//maxmum lenght of input strings 
	int[] M; //array to keep individual input strings' lenghts
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

	public int retLen; 
	long retTime;	
	partSol Global_finalSol ;

	class partSol// Solution
	{
		int[] X;
		int[] P;//the respected positions
		double h_pr;//for the heuristics based on probability 

		partSol(int pN, int pMaxLen)
		{
			X= new int[2*pMaxLen];
			P=new int[pN];
		}
	}//partSol


	//constructor
	IBS_SCS    (int Pn,  String pDataFile,  String string_type) throws IOException

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
			System.out.println("ERROOOOOOOR: real sigSize=" + Alphabet.length() + " as opposed to " + sigSize + " in " + pDataFile);
			return;
		}
		n=i+1;
		f.close();  
	}//constructor
	//-------------------------------------------------

	//-------------------------------------------------
	private void BShPr(int pBeamSize,int NUM_BEST, partSol pFinalSol)
	{
		int i2,j,k,kx,k2,t;
		partSol[] B=new partSol[pBeamSize];//the current list of (partial) solutions
		partSol[] C=new partSol[pBeamSize*sigSize];//for the children
		int nB,nC;//current size of the lists
		boolean[] flgCandids=new boolean[pBeamSize*sigSize];//for the children
		int curLen;
		boolean tmpChildUsable,flgAnswer;

		int tmpMaxHInx,tmpCount;

		double tmpLog2= Math.log(2);
		int MAX_TARGET= (int) (1+(m *Math.log(sigSize))/tmpLog2);


		double[][] Pr= new double[MAX_TARGET][m+1];

		double tmpCoeff=1;
		int jm,jk,tmpTarget;

		int minR, maxR, cut; 

		int  ib, tmpIb; int selectedChindIndex  ;
		int[] bestInc= new int[NUM_BEST];
		int[] tmpDFCounts= new int[NUM_BEST];
		boolean flgDF;
		for(ib=0;ib<NUM_BEST;ib++)
			tmpDFCounts[ib]=0;

		for(jm=0;jm<MAX_TARGET;jm++)
		{
			for(jk=0;jk<m+1;jk++)
			{
				if(jk==0)
				{
					Pr[jm][jk]= tmpCoeff;//1;
				}
				else
				{
					if(jk>jm)
					{
						Pr[jm][jk]= 0;
					}
					else
					{
						Pr[jm][jk]= tmpCoeff * 
						((1/(double)sigSize)*Pr[jm-1][jk-1] + ((sigSize-1)/(double)sigSize)*Pr[jm-1][jk]);
					}
				}
			}
		}

		for(k=0;k<B.length;k++)
			B[k]=new partSol(n,sigSize*m);
		for(k=0;k<C.length;k++)
			C[k]=new partSol(n,sigSize*m);
		for(i2=0;i2<n;i2++)
			B[0].P[i2]=-1;
		nB=1;
		nC=0;
		curLen=0;
		while(true) //until no more expansion possible
		{
			minR= m+1;//for initilization
			maxR= -1;//for initilization
			for(k=0;k<nB;k++)//for each candidate partial solution
			{
				for(t=0;t<sigSize;t++)//for each letter in Sigma
				{
					//check usability & answership
					tmpChildUsable=false;
					flgAnswer=true;
					for(i2=0;i2<n;i2++)
					{
						if( B[k].P[i2]<M[i2]-1 && S[i2].s[B[k].P[i2]+1]==t )
						{
							C[nC].P[i2]= B[k].P[i2]+1;
							tmpChildUsable=true;
						}
						else
						{
							C[nC].P[i2]= B[k].P[i2];
						}
						if(C[nC].P[i2] < M[i2]-1)
							flgAnswer=false;
						if(M[i2]-C[nC].P[i2]-1 < minR) 
							minR= M[i2]-C[nC].P[i2]-1;
						if(M[i2]-C[nC].P[i2]-1 > maxR)
							maxR= M[i2]-C[nC].P[i2]-1;
					}
					if(tmpChildUsable)// NOTE the order
					{
						for(kx=0;kx<curLen;kx++)
							C[nC].X[kx]= B[k].X[kx];
						C[nC].X[curLen]= t;//NOTE: starting from 0
						nC++;
					}
					if(flgAnswer)// NOTE the order
					{
						//prepare the final answer
						retLen=curLen+1;
						for(j=0;j<=curLen;j++)
							pFinalSol.X[j]=C[nC-1].X[j];
						for(i2=0;i2<n;i2++)
							pFinalSol.P[i2]=C[nC-1].P[i2];
						return;
					}
				}//for t
			}//for k

			if(nC>0)
			{
				tmpTarget= (int) ((maxR) * Math.log(sigSize)/tmpLog2);
				if(tmpTarget>=MAX_TARGET-1)	tmpTarget=  (maxR > MAX_TARGET-1)? maxR:(MAX_TARGET-1);
				if(tmpTarget<=maxR) tmpTarget=maxR;

				//determine cut
				cut=0;
				if(maxR > 100)
				{
					cut= maxR-100;
					if(minR-cut <1)
						cut= minR-1;
					tmpTarget= (int) ((maxR-cut) * Math.log(sigSize)/tmpLog2);
					if(tmpTarget>=MAX_TARGET-1)	tmpTarget=  ((maxR-cut) > MAX_TARGET-1)? (maxR-cut):(MAX_TARGET-1);
					if(tmpTarget <= (maxR-cut)) tmpTarget= (maxR-cut);
				}

				for(ib=0;ib<NUM_BEST;ib++)
					bestInc[ib]= 0;
				selectedChindIndex = 0;

				for(k2=0;k2<nC;k2++)
				{
					C[k2].h_pr=0;
					for(i2=0;i2<n;i2++)
						C[k2].h_pr= C[k2].h_pr + Math.log(Pr[tmpTarget][M[i2]-1-C[k2].P[i2]  -  cut ] );

					tmpIb=0;
					while(tmpIb < selectedChindIndex && tmpIb<NUM_BEST && C[k2].h_pr <= C[bestInc[tmpIb]].h_pr)
						tmpIb++;
					if(tmpIb<NUM_BEST)
					{
						for(ib=NUM_BEST-1; ib>tmpIb ;ib--)
							bestInc[ib]= bestInc[ib-1];
						bestInc[tmpIb]= k2;
						selectedChindIndex++;
					}
				}//for k2

				for(k2=0;k2<nC;k2++)
				{
					flgCandids[k2]=true;
					//check for DF with respect to the best candidate(s) found so far
					for(ib=0; flgCandids[k2] && ib<NUM_BEST; ib++)
					{
						if(k2 != bestInc[ib] && flgCandids[bestInc[ib]]) 
						{
							flgDF=true;
							for(i2=0; flgDF && i2<n; i2++)
							{
								if(C[k2].P[i2] > C[bestInc[ib]].P[i2])
									flgDF=false;
							}
							if(flgDF)
							{
								flgCandids[k2]=false;
								tmpDFCounts[ib]++;
							}
						}
					}
				}

				//select (at-most) pBeamSize ones
				nB=Math.min(nC,pBeamSize);
				tmpCount=0;
				while(tmpCount<nB)
				{
					//initialise tmpMaxHInx to the first available (candidate) solution
					tmpMaxHInx= 0;
					while(tmpMaxHInx<nC && !flgCandids[tmpMaxHInx])
						tmpMaxHInx++;

					//look for a better one, if any
					for(k2=tmpMaxHInx+1;k2<nC;k2++)
						if(flgCandids[k2] && C[k2].h_pr > C[tmpMaxHInx].h_pr)
							tmpMaxHInx=k2;
					if(tmpMaxHInx < nC)
					{
						//select:
						for(j=0;j<=curLen;j++)//Note =
							B[tmpCount].X[j]=C[tmpMaxHInx].X[j];
						for(i2=0;i2<n;i2++)
							B[tmpCount].P[i2]=C[tmpMaxHInx].P[i2];
						flgCandids[tmpMaxHInx]= false;//to avoid reseltion
						tmpCount++;
					}
					else
						break;
				}
				nB=tmpCount;
				nC=0;
				curLen++;
			}//if nC>0
			else
			{				
				String strMes ="Error!!!!!!!!!!!!!!: Unexpected exit point";
				System.out.println(strMes);
				return;
			}//else nC>0
		}//while true
	}//BShPr

	//-------------------------------------------------
	public void mainMethod(int beam_size, int NUM_BEST)
	{
		int i,j,k;
		partSol finalSol= new partSol(n,sigSize*m);
		BShPr(beam_size,NUM_BEST,finalSol);
		Global_finalSol = finalSol;
		//test the answer:
		boolean flgOK=true;
		for(i=0; i<n && flgOK; i++)
		{
			k=0;
			for(j=0; j<M[i] && flgOK; j++)
			{
				while(k<retLen && S[i].s[j]!=finalSol.X[k])
					k++;
				if(k>=retLen)
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

}//IBS_SCS

