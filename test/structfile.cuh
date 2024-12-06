#ifndef STRUCTFILE_H_
#define STRUCTFILE_H_
#include <iostream>
#include <string>
using namespace std;
namespace STRUCTFILE
{
    struct FilePath
    {
        string DEMPath;
        string inR0Path;
        string inR1Path;
        string inR2Path;
        string inR3Path;
        string inR4Path;
        string inR5Path;
        string inH_initialPath;
        string inU_initialPath;
        string inV_initialPath;
        string inSoil_depthPath;
        string inmanningPath;
        string inINF_KsPath;
        string inINF_UPath;
        string inINF_OsPath;
        string inINF_OiPath;
        string inlandusePath;
        string outHPath;
        string outUPath;
        string outVPath;
        string outCPath;
        string outh_maxPath;
        // string outhu_maxPath;
        // string outhv_maxPath;
        string outhvel_maxPath;
        string outINF_TotalPath;
        string outPltPath;
        FilePath()
        {
            DEMPath = "./input/z.txt";
            inR0Path = "./input/channel.txt";
            inR1Path = "./input/R/R1.txt";
            inR2Path = "./input/R/R2.txt";
            inR3Path = "./input/R/R3.txt";
            inR4Path = "./input/R/R4.txt";
            inR5Path = "./input/R/R5.txt";
            inH_initialPath = "./input/Iniital/H_initial.txt";
            inU_initialPath = "./input/Iniital/U_initial.txt";
            inV_initialPath = "./input/Iniital/V_initial.txt";
            inmanningPath = "./input/manning.txt";
            inSoil_depthPath = "./input/soil_depth.txt";
            inSoil_depthPath = "./input/soil_depth.txt";
            inINF_KsPath = "./input/Landuse/INF_Ks.txt";
            inINF_UPath = "./input/Landuse/INF_U.txt";
            inINF_OsPath = "./input/Landuse/INF_Os.txt";
            inINF_OiPath = "./input/Landuse/INF_Oi.txt";
            inlandusePath = "./input/landuse.txt";
            outHPath = "./output/H";
            outUPath = "./output/U";
            outVPath = "./output/V";
            outCPath = "./output/C";
            outh_maxPath = "./output/h_max";
            // outhu_maxPath = "./output/hu_max";
            // outhv_maxPath = "./output/hv_max";
            outhvel_maxPath = "./output/hvel_max";
            outINF_TotalPath = "./output/INF_Total/INF_Total";
            outPltPath = "./output/TEST.plt";

            // DEMPath = "z.txt";
            // inRPath = "input\\R\\R.txt";
            // inR1Path = "input\\R\\R1.txt";
            // inR2Path = "input\\R\\R2.txt";
            // inR3Path = "input\\R\\R3.txt";
            // inR4Path = "input\\R\\R4.txt";
            // inR5Path = "input\\R\\R5.txt";
            // inH_initialPath = "input\\Iniital\\H_initial.txt";
            // inU_initialPath = "input\\Iniital\\U_initial.txt";
            // inV_initialPath = "input\\Iniital\\V_initial.txt";
            // inSoil_depthPath = "input\\Landuse\\Soil_depth.txt";
            // inINF_KsPath = "input\\Landuse\\INF_Ks.txt";
            // inINF_UPath = "input\\Landuse\\INF_U.txt";
            // inINF_OsPath = "input\\Landuse\\INF_Os.txt";
            // inINF_OiPath = "input\\Landuse\\INF_Oi.txt";
            // outHPath = "output\\H\\H";
            // outUPath = "output\\U\\U";
            // outVPath = "output\\V\\V";
            // outCPath = "output\\C\\C";
            // outh_maxPath = "output\\h_max";
            // outhu_maxPath = "output\\hu_max";
            // outhv_maxPath = "output\\hv_max";
            // outhvel_maxPath = "output\\hvel_max";
            // outINF_TotalPath = "output\\INF_Total\\INF_Total";
            // outPltPath = "output\\TEST.plt";
        }
    };

} // namespace STRUCTFILE
#endif