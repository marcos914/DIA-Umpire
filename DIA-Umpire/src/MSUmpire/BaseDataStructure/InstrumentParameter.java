/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.BaseDataStructure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class InstrumentParameter implements Serializable{
    private static final long serialVersionUID = 7563887811638875862L;

    public int Resolution;
    public float MS1PPM;
    public float MS2PPM;
    public float SNThreshold;
    public float MinMSIntensity;
    public float MinMSMSIntensity;
    public int NoPeakPerMin = 150;
    public float MinRTRange;
    public int StartCharge = 2;
    public int EndCharge = 5;
    public int MS2StartCharge = 2;
    public int MS2EndCharge = 4;
    public float MaxCurveRTRange=2f; 
    public float RTtol;
    public float MS2SNThreshold;
    public InstrumentType InsType;
    public int MaxNoPeakCluster=4;
    public int MinNoPeakCluster=2;
    public int MaxMS2NoPeakCluster=3;
    public int MinMS2NoPeakCluster=2;
    public boolean Denoise = true;
    public boolean EstimateBG = false;
    public boolean DetermineBGByID = false;
    public boolean RemoveGroupedPeaks = true;
    public boolean Deisotoping = false;
    public transient boolean BoostComplementaryIon=true; 
    public transient boolean AdjustFragIntensity=true;
    public int PrecursorRank = 25;
    public int FragmentRank = 300;
    public float RTOverlapThreshold = 0.1f;
    public float CorrThreshold = 0.1f;
    public float ApexDelta = 0.6f;
    public float SymThreshold = 0.3f;
    public int NoMissedScan = 1;
    public int MinPeakPerPeakCurve = 1;
    public float MinMZ=200;
    public int MinFrag=10;
    public transient float MiniOverlapP =  0.2f;
    public transient boolean CheckMonoIsotopicApex=false;
    public transient boolean DetectByCWT=true;
    public transient boolean FillGapByBK=true;
    public transient float IsoCorrThreshold=0.2f;
    public transient float RemoveGroupedPeaksCorr=0.3f;
    public transient float RemoveGroupedPeaksRTOverlap=0.3f;
    public transient float HighCorrThreshold=0.7f;
    public transient int MinHighCorrCnt=10;
    public transient int TopNLocal=6;
    public transient int TopNLocalRange=100;
    public transient float IsoPattern = 0.3f;
    public transient float startRT=0f;
    public transient float endRT=9999f;
    public transient boolean TargetIDOnly=false;
    public transient boolean MassDefectFilter=true;
    public transient float MinPrecursorMass=600f;
    public transient float MaxPrecursorMass=15000f;
    public transient boolean UseOldVersion=false;
    public transient float RT_window_Targeted=-1f;
    public transient int SmoothFactor = 5;
    public transient boolean DetectSameChargePairOnly=false;
    public transient float MassDefectOffset=0.1f;
    public transient int MS2PairTopN=5;    
    public transient boolean MS2Pairing=true;
    
    
    public void WriteParamSerialization(String mzXMLFileName) {
        try {
            Logger.getRootLogger().info("Writing parameter to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_params.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_params.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public static InstrumentParameter ReadParametersSerialization(String filepath) {
        if(! new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_params.ser").exists()){
            return null;
        }
        try {
            Logger.getRootLogger().info("Reading parameters from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_params.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_params.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            InstrumentParameter params = (InstrumentParameter) in.readObject();
            in.close();
            fileIn.close();
            return params;

        } catch (Exception ex) {
             Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
    }
 
    public InstrumentParameter(InstrumentType type) {
        InsType = type;
        SetParameter(type);
    }

    public static float CalcPPM(float valueA, float valueB) {
        return Math.abs(valueA - valueB) * 1000000 / valueB;
    }
    
    /**
     *
     * @param valueA
     * @param charge
     * @param ppm
     * @return
     */
    public static float GetMzByPPM(float valueA, int charge, float ppm) {
        float mwA = valueA * charge - charge * 1.00727f;
        float premass = mwA - (ppm * mwA / 1000000);
        return (premass + charge * 1.00727f) / charge;
    }

    public static float CalcSignedPPM(float valueA, float valueB) {
        return (valueA - valueB) * 1000000 / valueB;
    }

    /**
     *
     */
    public enum InstrumentType {
        Fusion,
        Q_TOF,
        TOF5600,
        Orbitrap,
        QExactive,
        MALDI,
    };

    private void SetParameter(InstrumentType type) {
        switch (type) {
            case Orbitrap: {
                MS1PPM = 10;
                MS2PPM = 20;
                SNThreshold = 3f;
                MinMSIntensity = 50f;
                MinMSMSIntensity = 20f;
                MinRTRange = 0.1f;
                MaxNoPeakCluster = 4;
                MinNoPeakCluster = 2;
                MaxMS2NoPeakCluster = 3;
                MinMS2NoPeakCluster = 2;
                MaxCurveRTRange = 2;
                RTtol = 0.1f;
                Resolution = 60000;
                break;
            }
            case QExactive: {
                MS1PPM = 10;
                MS2PPM = 20;
                SNThreshold = 3f;
                MS2SNThreshold = 3f;
                MinMSIntensity = 500f;
                MinMSMSIntensity = 100f;
                MinRTRange = 0.1f;
                MaxNoPeakCluster = 4;
                MinNoPeakCluster = 2;
                MaxMS2NoPeakCluster = 3;
                MinMS2NoPeakCluster = 2;
                MaxCurveRTRange = 1;
                RTtol = 0.1f;
                NoPeakPerMin = 150;
                Resolution = 30000;
                NoMissedScan = 2;
                Denoise = true;
                EstimateBG = false;
                RemoveGroupedPeaks = true;
                break;
            }
             case Fusion: {
                MS1PPM = 10;
                MS2PPM = 20;
                SNThreshold = 2f;
                MS2SNThreshold = 2f;
                MinMSIntensity = 100f;
                MinMSMSIntensity = 100f;
                MinRTRange = 0.1f;
                MaxNoPeakCluster = 4;
                MinNoPeakCluster = 2;
                MaxMS2NoPeakCluster = 3;
                MinMS2NoPeakCluster = 2;
                MaxCurveRTRange = 2f;
                RTtol = 0.1f;
                NoPeakPerMin = 150;
                Resolution = 50000;
                NoMissedScan = 1;
                Denoise = true;
                EstimateBG = false;
                RemoveGroupedPeaks = true;
                break;
            }
            case Q_TOF: {
                MS1PPM = 40;
                MS2PPM = 100;
                SNThreshold = 3f;
                MinMSIntensity = 10f;
                MinMSMSIntensity = 10f;
                MinRTRange = 0.3f;
                MaxNoPeakCluster = 5;
                MinNoPeakCluster = 3;
                MaxMS2NoPeakCluster = 3;
                MinMS2NoPeakCluster = 2;
                MaxCurveRTRange = 2;
                Resolution = 17000;
                RTtol = 0.1f;
                break;
            }
            case TOF5600: {
                MS1PPM = 30;
                MS2PPM = 40;
                SNThreshold = 2f;
                MS2SNThreshold = 2f;
                MinMSIntensity = 5f;
                MinMSMSIntensity = 1f;
                MinRTRange = 0.1f;
                MaxNoPeakCluster = 4;
                MinNoPeakCluster = 2;
                MaxMS2NoPeakCluster = 4;
                MinMS2NoPeakCluster = 2;
                MaxCurveRTRange = 1.5f;
                Resolution = 17000;
                RTtol = 0.1f;
                Denoise = true;
                EstimateBG = true;
                RemoveGroupedPeaks = true;
                break;
            }
        }
    }
}
