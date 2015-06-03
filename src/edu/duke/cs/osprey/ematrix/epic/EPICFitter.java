/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix.epic;

/**
 *
 * This object calculates an EPIC fit for an objective function
 * over continuously flexible degrees of freedom
 * 
 * @author mhall44
 */
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;


//this class generates EPIC polynomial fits for a given RC (intra+shell energy) or pair (pairwise energy)
//it can try various orders until it gets a good one
public class EPICFitter {

    EPICSettings es;
    
    //template for the generated EPoly objects
    public int numDOFs;
    DoubleMatrix1D DOFmax;
    DoubleMatrix1D DOFmin;
    DoubleMatrix1D center;
    double minE;
    
    public EPoly PCTemplate = null;//template to use for EPolyPCs
    //needs to be set (by quadratic fitting) before trying to make any of these
    
    MolecEObjFunction objFcn;//the energy as a function of DOF values...this
    //is what we're trying to approximate
    
    static int sampPerParam = 10;
    
    public EPICFitter ( MolecEObjFunction mof, EPICSettings eset,
            DoubleMatrix1D cen, double me ) {
        //given the CCDMinimizer used to minimize for a rotamer pair (or intra+shell)
        //construct an EpicFitter for it
        
        objFcn = mof;
        
        numDOFs = mof.getNumDOFs();
        DOFmax = mof.getConstraints()[1];
        DOFmin = mof.getConstraints()[0];

        es = eset;
        center = cen;
        minE = me;
    }

    
    
    
    
    public EPoly doFit(FitParams fp){

        int numParams = SeriesFitter.getNumParams(numDOFs, false, fp.order);
        int numSamples = sampPerParam * numParams;
        
        EPolyPC PCFit = null;


        if(fp.PCOrder>fp.order){//need to do PC fit

            PCFit = new EPolyPC(PCTemplate,fp.order,fp.PCOrder,fp.PCFac);

            int numPCs = SeriesFitter.countTrue(PCFit.isPC);

            fp.numPCParams = 0;
            
            for(int n=fp.order+1; n<=fp.PCOrder; n++)
                fp.numPCParams += SeriesFitter.getNumParamsForOrder(numPCs,n);
            
            numParams += fp.numPCParams;

            numSamples = sampPerParam*numParams;
        }


        
        //PARAMETER CAP: save time
        if(numParams>2000){
            System.out.println("ABORTING EPICFITTER.DOFIT BECAUSE THERE ARE TOO MANY PARAMETERS: "+numParams);
            return null;
        }

        
        DoubleMatrix1D[] sampRel = new DoubleMatrix1D[numSamples];
        DoubleMatrix1D[] sampAbs = new DoubleMatrix1D[numSamples];
        double trueVal[] = new double[numSamples];

        generateSamples(numSamples,sampRel,sampAbs,trueVal,numSamples/2);

        //no point doing high1S if all samples are below bCutoff
        boolean allBelowCutoff = true;
        for(int s=0; s<numSamples; s++){
            if(trueVal[s]>=es.EPICThresh1){
                allBelowCutoff = false;
                break;
            }
        }
        
        SAPE sapeTerm = null;
        double bCutoffs[] = new double[numSamples];
        double bCutoffs2[] = new double[numSamples];
        double baseShift = 0;//SVE at center
        
        if(fp.SAPECutoff==0){
            Arrays.fill(bCutoffs, es.EPICThresh1);
            Arrays.fill(bCutoffs2, es.EPICThresh2);
        }
        //if using double bCutoff for provability also set bCutoffs array
        //according to that here!
        else{//>0            
            sapeTerm = new SAPE(objFcn,fp.SAPECutoff,sampAbs);
            
            baseShift = sapeTerm.getEnergyStandalone(center);//SAPE contribution at center
            
            for(int s=0; s<numSamples; s++){
                double shift = sapeTerm.getEnergyStandalone(sampAbs[s]);
                trueVal[s] -= shift - baseShift;
                bCutoffs[s] = es.EPICThresh1 - shift + baseShift;//this keeps bCutoff at the same place
                //accounting for the sparse VDW energy
                
                bCutoffs2[s] = es.EPICThresh2 - shift + baseShift;
            }
        }
        

        double weights[] = null;
            //trying out weighted least squares
        weights = new double[numSamples];
        for(int s=0; s<numSamples; s++){
            if(trueVal[s]>1)
                weights[s] = 1/trueVal[s];
            else
                weights[s] = 1;
        }


        double lambda = 0;

        double[] seriesCoeffs = null;//used for order>2
        EPoly ans = null;

        if(fp.PCOrder>fp.order){//use principal components
            //fit in PC basis
            DoubleMatrix1D ySamp[] = new DoubleMatrix1D[numSamples];
            for(int s=0; s<numSamples; s++)
                ySamp[s] = PCFit.toPCBasis(sampRel[s]);

            if(allBelowCutoff){
                System.out.println("Analytical:");
                PCFit.coeffs = SeriesFitter.fitSeries(ySamp,trueVal,weights,lambda,
                    false,PCFit.fullOrder,PCFit.PCOrder,PCFit.isPC,false,null,null);
            }
            else {
                //ITERATIVE
                PCFit.coeffs = SeriesFitter.fitSeriesIterative(ySamp,trueVal,weights,lambda,
                    false,PCFit.fullOrder,bCutoffs,bCutoffs2,PCFit.PCOrder,PCFit.isPC);
            }
            
            ans = PCFit;
        }
        else{
            
            if(allBelowCutoff){
                System.out.println("Analytical:");
                seriesCoeffs = SeriesFitter.fitSeries(sampRel,trueVal,weights,lambda,
                        false,fp.order);
            }
            else {
                seriesCoeffs = SeriesFitter.fitSeriesIterative(sampRel,trueVal,weights,
                        lambda,false,fp.order,bCutoffs,bCutoffs2,fp.order,null);
            }
            
            ans = new EPoly(numDOFs, objFcn.getDOFs(), DOFmax, DOFmin, center,
                    minE, seriesCoeffs, fp.order);
        }
        
        
        //add SAPE term
        if(fp.SAPECutoff>0){
            ans.sapeTerm = sapeTerm;
            ans.baseSAPE = baseShift;
        }
        
        
        
        //DEBUG!!!!!
            //checking training set mean residual!
        //similar to crossValidateSeries
        /*if(fp.PCOrder<=fp.order){
            double meanResidual = 0;
            double weightSum = 0;


            for(int s=0; s<numSamples; s++){

                DoubleMatrix1D x = sampAbs[s];

                double realVal = ccdMin.objFcn.getValue(x) - ans.minE;
                //actual energy, relative to voxel minimum

                double sampBCutoff = es.EPICThresh1;
                double sampBCutoff2 = es.EPICThresh2;
                double serVal = ans.evaluate(x,false);

                if(ans.sve!=null){
                    ans.sveOF.setDOFs(x);
                    double shift = ans.sve.getEnergy();
                    realVal -= shift - baseShift;
                    sampBCutoff -= shift - baseShift;
                    sampBCutoff2 -= shift - baseShift;
                    serVal -= shift - baseShift;
                    //we subtract off the SVE contribution (shift-baseShift) from everything
                    //(-baseShift because we set the SVE zero point at center)
                }

                double weight = 1;
                if(realVal>1)
                    weight = 1/realVal;//1/Math.min(realVal,sampBCutoff);
                weightSum += weight;


                if(realVal>=sampBCutoff){
                    if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, false)){
                        meanResidual += weight * (serVal-sampBCutoff)*(serVal-sampBCutoff);
                    }
                    if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true)){
                        meanResidual += weight * (realVal-serVal)*(realVal-serVal);
                    }
                }
                else
                    meanResidual += weight*(realVal-serVal)*(realVal-serVal);
            }

            meanResidual /= weightSum;//numSamples;
            System.out.println("CHECK TRAINING SET MEAN RESIDUAL:"+meanResidual);
        }
        */    //DEBUG!!!!
        
        ans.fitDescription = fp.getDescription();
        
        return ans;
    }



    public double crossValidateSeries(EPoly fit, FitParams fp){
        //(note this is much like generateSamples!!)
        //return mean residual
        //used in addFitSeries and augmentSeries
        
        if(fit==null)//parameter cap will do this
            return Double.POSITIVE_INFINITY;//obviously an absent fit is no good
        
        double meanResidual = 0;
        double weightSum = 0;

        double baseShift = 0;
        if(fit.sapeTerm!=null){
            baseShift = fit.sapeTerm.getEnergyStandalone(center);
        }
        
        
        /*double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = DOFmax.get(dof) - center.get(dof);
            relMin[dof] = DOFmin.get(dof) - center.get(dof);
        }*/
        
        int numSamples = sampPerParam * fp.numParams();
        
        double trueVal[] = new double[numSamples];
        DoubleMatrix1D sampRel[] = new DoubleMatrix1D[numSamples];
        DoubleMatrix1D sampAbs[] = new DoubleMatrix1D[numSamples];
        generateSamples(numSamples,sampRel,sampAbs,trueVal,numSamples/2);
        
        for(int s=0; s<numSamples; s++){

            /*DoubleMatrix1D dx = DoubleFactory1D.dense.make(numDOFs);
            DoubleMatrix1D x = DoubleFactory1D.dense.make(numDOFs);
            for(int dof=0; dof<numDOFs; dof++){
                double top = relMax[dof];
                double bottom = relMin[dof];
                dx.set(dof, bottom + Math.random()*(top-bottom));
                x.set(dof, center.get(dof)+dx.get(dof));
            }

            double realVal = ccdMin.objFcn.getValue(x) - fit.minE;
            //actual energy, relative to voxel minimum
            */
            DoubleMatrix1D x = sampAbs[s];
            double realVal = trueVal[s];

            double sampBCutoff = es.EPICThresh1;
            double sampBCutoff2 = es.EPICThresh2;
            double serVal = fit.evaluate(x,false,false);
            
            if(fit.sapeTerm!=null){
                double shift = fit.sapeTerm.getEnergyStandalone(x);
                realVal -= shift - baseShift;
                sampBCutoff -= shift - baseShift;
                sampBCutoff2 -= shift - baseShift;
                serVal -= shift - baseShift;
                //we subtract off the SVE contribution (shift-baseShift) from everything
                //(-baseShift because we set the SVE zero point at center)
            }
            
            double weight = 1;
            if(realVal>1)
                weight = 1/realVal;//1/Math.min(realVal,sampBCutoff);
            weightSum += weight;

            
            /*if(serVal<sampBCutoff||realVal<sampBCutoff)
                //meanResidual += (realVal-bv)*(realVal-bv);
                meanResidual += weight * (realVal-serVal)*(realVal-serVal);
            else if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true))//bCutoff2 upper restraint
                meanResidual += weight*(realVal-serVal)*(realVal-serVal);*/
            if(realVal>=sampBCutoff){
                if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, false)){
                    meanResidual += weight * (serVal-sampBCutoff)*(serVal-sampBCutoff);
                }
                if(SeriesFitter.isRestraintTypeActive(realVal, serVal, sampBCutoff, sampBCutoff2, true)){
                    meanResidual += weight * (realVal-serVal)*(realVal-serVal);
                }
            }
            else
                meanResidual += weight*(realVal-serVal)*(realVal-serVal);
        }

        meanResidual /= weightSum;//numSamples;
        System.out.println("CV MEAN RESIDUAL:"+meanResidual);
        
        //Let's return the mean residual
        return meanResidual;
        //END CV
    }



    static void analyzeLSBRecord(ArrayList<double[]> LSBRecord){
        //Given the triples (old lower bound, LSB lower bound, minimized energy)
        //for a bunch of conformations
        //provide statistics on them


        double binMaxs[] = {0.1,0.5,0.75,0.9,0.99,1,1.01,1.1,1.5,Double.POSITIVE_INFINITY};
        int numBins = binMaxs.length;
        double binMins[] = new double[numBins];
        binMins[0] = Double.NEGATIVE_INFINITY;
        System.arraycopy(binMaxs,0,binMins,1,numBins-1);


        double avgSR = 0;//average of slack recovery fractions
        int binCounts[] = new int[numBins];//how many conformations are in each bin
        int numOverHundredth = 0;//Number of LSB bounds > 0.01 kcal/mol over true energy
        int numOverTenth = 0;//>0.1 kcal/mol over
        int numOverHalf = 0;//>0.5 kcal/mol over
        int numEnum=0;//how many conformations we'd need to enumerate
        //using the LSB lower bound



        
        double minMinE = Double.POSITIVE_INFINITY;//Actual GMEC energy
        for(double[] rec: LSBRecord)
            minMinE = Math.min(minMinE,rec[2]);


        double avgTimeRat = 0;

        for(double[] rec : LSBRecord){

            double slackRecovered = (rec[1]-rec[0])/(rec[2]-rec[0]);
            //fraction of slack recovered by LSB

            avgSR += slackRecovered;

            if(rec[1]<=minMinE)
                numEnum++;
            if(rec[1]>rec[2]+0.01)
                numOverHundredth++;
            if(rec[1]>rec[2]+0.1)
                numOverTenth++;
            if(rec[1]>rec[2]+0.5)
                numOverHalf++;

            for(int bin=0; bin<numBins; bin++){
                if(slackRecovered<=binMaxs[bin]&&slackRecovered>binMins[bin]){
                    binCounts[bin]++;
                    break;
                }
            }
            
            avgTimeRat += rec[3]/rec[4];
        }

        avgSR /= LSBRecord.size();
        avgTimeRat /= LSBRecord.size();
        
        System.out.println("ANALYSIS OF LSB:");
        System.out.println("Total conformation count: "+LSBRecord.size());
        System.out.println("Average minimization time ratio (normal/EPIC): "+avgTimeRat);
        System.out.println("Average slack recovery fraction: "+avgSR);
        System.out.println(numOverHundredth+" LSBs > 0.01 over true E; "+numOverTenth+" >0.1 over, "
                +numOverHalf+" >0.5 over");
        //System.out.println(numEnum+" need to be enumerated based on LSBs");
        System.out.println("Bin_max Bin_count");
        for(int bin=0; bin<binMaxs.length; bin++){
            System.out.println(binMaxs[bin]+" "+binCounts[bin]);
        }
    }
    
    



    void generateSamples(int numSamples,
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal, int maxOverCutoff){

        //Generate samples relative to startVec (sampRel) and absolute (sampAbs)
        //and give their energies (trueVal), relative to baseE
        //we'll allow at most maxOverCutoff samples above es.EPICThresh1
        //to prevent overfitting

        //by default, uniform voxel sampling is used
        //shift to Metropolis if we end up drawing numSamples/2 samples over bCutoff
        //and by that time less than maxUniFrac*numSamples/2 have been drawn under bCutoff
        //double maxUniFrac = 0.05;
        
        //Let's instead say if uniform draw fails enough times (maxFailCount) in a row we bail out and do Metropolis
        int maxFailCount = 10000;//big because Metropolis is expensive
        
        //DEBUG!!  To test Metropolis
        //generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,0);
        //return;
        //DEBUG!!!
        
        
        
        int countOverCutoff = 0;//how many of our samples are over the cutoff
                        
        double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = DOFmax.get(dof) - center.get(dof);
            relMin[dof] = DOFmin.get(dof) - center.get(dof);
        }


        
        for(int s=0; s<numSamples; s++){

            if(countOverCutoff<maxOverCutoff){//normal draw
                uniformVoxelSample(s,sampRel,sampAbs,trueVal,objFcn,relMin,relMax);
                if(trueVal[s]>es.EPICThresh1)
                    countOverCutoff++;
            }
            else {//force sub-threshold draw
                /*if(s-countOverCutoff<maxUniFrac*countOverCutoff){//not getting enough sub-threshold samples
                    //switch to Metropolis
                    generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,s);
                    break;
                }
                else {//getting enough some-threshold samples to keep up the uniform method
                    */
                    int failCount = 0;//how many uniform samples ended up above the cutoff in a rows
                    
                    do {
                        uniformVoxelSample(s,sampRel,sampAbs,trueVal,objFcn,relMin,relMax);
                        
                        if(failCount>maxFailCount){//bail and do Metropolis
                            generateSamplesMetropolis(numSamples,sampRel,sampAbs,trueVal,s);
                            System.out.println("Drew "+numSamples+" samples of which "+countOverCutoff+" are over bCutoff");
                            return;
                        }
                        failCount++;
                            
                    } while(trueVal[s]>es.EPICThresh1);
                //}
            }
        }
        
        System.out.println("Drew "+numSamples+" samples of which "+countOverCutoff+" are over bCutoff");
    }
            

    void generateSamplesMetropolis(int numSamples,
            DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal, int startSampNum){
        //Use Metropolis (in SubThreshSampler) to get samples
        //way less efficient than uniform sampling if much of the voxel is below EPICThresh1
        //but can hopefully deal with very small hypervolumes being below EPICThresh1
        //this function can start drawing at 
        //so we can have uniformly drawn samples before it (e.g. the above-the-cutoff ones)
        
        SubThreshSampler sts = new SubThreshSampler(es.EPICThresh1+minE,objFcn,DOFmin,DOFmax);
        
        sts.burnIn(center);//start at center because it's known to be below the threshold
        
        for(int s=startSampNum; s<numSamples; s++){
            sampAbs[s] = sts.nextSample();
            trueVal[s] = objFcn.getValue(sampAbs[s]) - minE;
            sampRel[s] = sampAbs[s].copy();
            sampRel[s].assign(center,Functions.minus);
        }
    }
    
    void uniformVoxelSample(int s, DoubleMatrix1D[] sampRel, DoubleMatrix1D[] sampAbs, double[] trueVal,
            ObjectiveFunction of, double[] relMin, double[] relMax){
        //Draw sample # s for generateSamples, filling in trueVal, sampRel, and sampAbs
        //Generate vector relative to minimum
        DoubleMatrix1D dx = DoubleFactory1D.dense.make(numDOFs);
        //and absolute
        DoubleMatrix1D x = DoubleFactory1D.dense.make(numDOFs);

        for(int dof=0; dof<numDOFs; dof++){
            double top = relMax[dof];
            double bottom = relMin[dof];

            dx.set(dof, bottom + Math.random()*(top-bottom));
            x.set(dof, center.get(dof)+dx.get(dof));
        }

        trueVal[s] = of.getValue(x) - minE;

        sampRel[s] = dx;
        sampAbs[s] = x;
    }
    
    
    public FitParams raiseFitOrder(FitParams fp){
        //given a FitParams object, return a higher-order one
        //to follow the standard EPIC progression of fit orders
        //(modified by useSVE, usePC as needed)
        //all these fits are meant to be performed iteratively (with thresholds) though
        //analytical evaluation will kick in if no samples exceed the first threshold
        //this function is intended to be used repetitively until a good enough fit is found, starting
        //with a FitParams.quadratic()
        
        //since these are standard EPIC fits we don't includeConst
        
        
        if(es.useSAPE){
            //we try SVE fits after non-PC, polynomial-only fits of the same order
            if(fp.SAPECutoff==0 && fp.order>=fp.PCOrder){
                double cutoff = 3;
                if(fp.order>2)//4 or 6
                    cutoff = 4;
                return new FitParams(numDOFs,fp.order,0,fp.order,false,cutoff);
            }
        }
                
        
        if(es.usePC && fp.order<6){
            //after quadratic or quartic fits without principal components (SVE or not)
            //we try first PCFac=0.1 and then 0.01
            if(fp.PCOrder==fp.order)//start PC with 0.1
                return new FitParams(numDOFs,fp.order,0.1,fp.order+2,false,0);
            else if(fp.PCFac==0.1)//move on to 0.01
                return new FitParams(numDOFs,fp.order,0.01,fp.order+2,false,0);
        }
        
        
        //if we get here we need to raise the main degree of the polynomial
        //we go 2, 4, 6, and then give up
        
        if(fp.order>=6)//give up
            return null;
        else//raise order by 2
            return new FitParams(numDOFs,fp.order+2,0,fp.order+2,false,0);
    }
    
    
    
    public EPoly blank(){
        //return a EPoly on no continuous degrees of freedom
        EPoly ans = new EPoly(numDOFs, objFcn.getDOFs(), DOFmax, DOFmin, center,
                minE, null, 2);
        //arbitrarily calling it quadratic (doesn't matter since no variables in polynomial)
        ans.fitDescription = "No DOFs";
        return ans;
    }

    
    
    void makeVoxelFigureData(EPoly ep){
        //This is a manually configured function that can be called from RotamerSearch.compEPICFit once fitting is completed
        System.out.println("Making voxel figure data!");
        
        
        System.out.println("minE: "+ep.minE);
        DoubleMatrix1D x = ep.center.copy();
        
        
        //this is for 1-D
        //if called on a >1-D voxel will stay at center wrt other dims
        if(ep.numDOFs!=2){
            System.out.println("chi trueval fitval");

            for(double chi=ep.DOFmin.get(0); chi<=ep.DOFmax.get(0); chi++){
                x.set(0, chi);
                double trueVal = objFcn.getValue(x);
                double fitVal = ep.evaluate(x, true, false);
                System.out.println(chi+" "+trueVal+" "+fitVal);
            }
        }
        else {
            //2-D version:
            System.out.println("chi1 chi2 trueval fitval");

            for(double chi1=ep.DOFmin.get(0); chi1<=ep.DOFmax.get(0); chi1++){
                for(double chi2=ep.DOFmin.get(1); chi2<=ep.DOFmax.get(1); chi2++){
                    x.set(0, chi1);
                    x.set(1, chi2);
                    double trueVal = objFcn.getValue(x);
                    double fitVal = ep.evaluate(x, true, false);
                    System.out.println(chi1+" "+chi2+" "+trueVal+" "+fitVal);
                }
            }
        }
            
        
        //We'll only want this done once, in a special OSPREY run
        System.out.println("Voxel figure data complete.");
        System.exit(0);
    }
    
}

