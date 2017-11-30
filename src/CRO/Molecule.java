package CRO;


/**
 * This class describes the molecule in RCCRO. An instance of this class represent
 * a molecule in the container. The elementary reactions that a molecule can be
 * involved is also implemented here.
 * @author C M Khaled Saifullah
 */

import java.util.Vector;

public class Molecule 
{
    /**
     * The total number of hits.
     */
    private int numHit;
    /**
     * The index of the last local optimum hit.
     */
    private int minHit;
    /**
     * The molecular structure this molecule holds.
     */
    private Vector<Integer> mols;
    /**
     * The potential and kinetic energy this molecule holds.
     */
    private double PE, KE;
    
	/**
     * The local minimum this molecule previously reached.
     */
    private double localMin;
    
    /**
     * The molecular structure this minimum potential energy.
     */
    private Vector<Integer> minMols;

	/**
     * This function is the class contructor.
     * @param mols Molecular Structure of the molecule.
     * @param initialKE inital KE .
     */
    public Molecule()
    {
    }
    
    /**
     * This function is a local update in the molecule-wide range.
     */
    public int update(int l) 
    {
    	//System.out.println(l);
        if(l==5)
        	l = 1;
        else if(l==10)
        	l = 5;
        else if(l==50)
        	l = 10;
        else if(l==100)
        	l = 10;
        else if(l==500)
        	l = 1000;
        else if(l==1000)
        	l = 1000;
        else if(l==5000)
        	l = 1000;
        
        return l;
    }
    
    public int update(int l,int status) 
    {
        if(l==100)
        	l = 3;
        else if(l==500)
        	l = 7;
        else if(l==1000)
        	l = 10;
        else if(l==5000)
        	l = 20;
        
        return l;
    }
    
    
    
    public int getNumHit() {
		return numHit;
	}
	public void setNumHit(int numHit) {
		this.numHit = numHit;
	}
	public int getMinHit() {
		return minHit;
	}
	public void setMinHit(int minHit) {
		this.minHit = minHit;
	}
	public Vector<Integer> getMols() {
		return mols;
	}
	public void setMols(Vector<Integer> mols) {
		this.mols = mols;
	}
	public double getPE() {
		return PE;
	}
	public void setPE(double pE) {
		PE = pE;
	}
	public double getKE() {
		return KE;
	}
	public void setKE(double kE) {
		KE = kE;
	}
	public double getLocalMin() {
		return localMin;
	}
	public void setLocalMin(double localMin) {
		this.localMin = localMin;
	}
	
    public Vector<Integer> getMinMols() {
		return minMols;
	}

	public void setMinMols(Vector<Integer> minMols) {
		this.minMols = minMols;
	}


}
