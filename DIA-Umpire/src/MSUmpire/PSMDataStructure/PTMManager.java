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
package MSUmpire.PSMDataStructure;

import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.io.SerializationUtils;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import org.apache.log4j.Logger;
import org.xmlpull.v1.XmlPullParserException;

/**
 * PTM library manager from compomics
 * @author Chih-Chiang Tsou
 */
public class PTMManager {

    private PTMFactory ptmFactory;
    private static PTMManager pTMManager;
    private String tempptmfile="ptmFactory-3.28.24.cus";
    private String tmpfilefolder=System.getProperty("user.home") + "/.compomics";

    public static PTMManager GetInstance() throws XmlPullParserException, IOException {
        if (pTMManager == null) {
            pTMManager = new PTMManager();
        }
        return pTMManager;
    }
       
    private PTMManager() throws XmlPullParserException, IOException { 
        //PTMManager.GetInstance();
        
        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/mods.xml");
        String tmpfile = "mods.xml";
        InputStreamToFile convert = new InputStreamToFile();

        File ptmFile = convert.GetFile(is, tmpfile);
        PTMFactory.setSerializationFolder(tmpfilefolder);
        
        PTMFactory.getInstance().clearFactory();
        PTMFactory.getInstance().importModifications(ptmFile, false, true);
        ptmFactory = PTMFactory.getInstance();   
        if (!ptmFactory.getDefaultModifications().isEmpty()) {
            SaveTempFile();
        } else {
            Logger.getRootLogger().error("Modification map file is empty");
        }
    }

    private void SaveTempFile() throws IOException {
        File factoryFile = new File(tmpfilefolder, tempptmfile);
        if (!factoryFile.getParentFile().exists()) {
            factoryFile.getParentFile().mkdir();
        }
        SerializationUtils.writeObject(ptmFactory, factoryFile);
    }
    
    
    
    public void ImportUserMod(String file) throws XmlPullParserException, IOException{
        File usermod=new File(file);        
        if (usermod.exists()) {
            ptmFactory.importModifications(usermod, true,false);
            if (!ptmFactory.getDefaultModifications().isEmpty()) {
                SaveTempFile();
            }
            else{
                Logger.getRootLogger().error("Modification map file is empty");;
            }
        }        
    }
    
    
    public static ArrayList<ModificationMatch> TranslateModificationString(String ModificationString) {        
        ArrayList<ModificationMatch> modlist = new ArrayList<>();
        final java.util.regex.Pattern p = java.util.regex.Pattern.compile("(.+?)\\((\\d+)\\)");
        if (ModificationString != null && !"".equals(ModificationString)) {
            for (final String mod : ModificationString.split(";")) {
                final java.util.regex.Matcher m = p.matcher(mod);
                m.lookingAt();
                final String ptmstring = m.group(1);
                final int site = Integer.parseInt(m.group(2));
                modlist.add(new ModificationMatch(ptmstring, true, site));
            }
        }
        return modlist;
    }
    
    
    public void AddPTM(PTM ptm) throws IOException{
        ptmFactory.addUserPTM(ptm);
        ptmFactory.saveFactory();
    }
    
    public PTM GetPTM(String AA, float massdiff) {

        double smallmassdiff = Double.MAX_VALUE;
        PTM smallestdiffptm = null;
        for (int i = 0; i < ptmFactory.getPTMs().size(); i++) {
            String name = ptmFactory.getPTMs().get(i);
            PTM ptm = ptmFactory.getPTM(name);
            boolean sitecorrect = false;
            if (("C-term".equals(AA) && name.toLowerCase().contains("c-term")) || ("N-term".equals(AA) && name.toLowerCase().contains("n-term"))) {
                sitecorrect = true;
            }
            if (ptm.getPattern() != null) {
                for (Character residue : ptm.getPattern().getAminoAcidsAtTarget()) {
                    if (String.valueOf(residue).equals(AA)) {
                        sitecorrect = true;
                    }
                }
            }
            if (sitecorrect) {
                double diff = Math.abs(ptm.getMass() - massdiff);
                if (diff < 0.5f) {
                    if (diff < smallmassdiff) {
                        smallmassdiff = diff;
                        smallestdiffptm = ptm;
                    }
                }
            }
        }
        return smallestdiffptm;
    }
    
    public PTM GetPTMByName(String modname) {

        for (int i = 0; i < ptmFactory.getPTMs().size(); i++) {
            String name = ptmFactory.getPTMs().get(i);
            PTM ptm = ptmFactory.getPTM(name);
            if (ptm.getName()== null ? modname == null : ptm.getName().equals(modname)) {
                return ptm;
            }
        }
        return null;
    }        
}
