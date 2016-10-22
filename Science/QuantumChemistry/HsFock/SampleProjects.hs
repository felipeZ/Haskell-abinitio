
{-|
Module: Science.QuantumChemistry.HsFock.SampleProjects
Description: Sample Data to Felipe
Copyright: @2013 Felipe Zapata, Angel Alvarez
           @2016 Felipe Zapata
The HaskellFock SCF Project
-}

module Science.QuantumChemistry.HsFock.SampleProjects where

import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.HartreeFock.HartreeFock
import Science.QuantumChemistry.HsFock.Project
import Science.QuantumChemistry.BasisSet.NormalizeBasis (normGlobal)


project :: String -> String -> ProjectData
project "HHe+" "STO-3G" =
    PD {
          pLabel           = "HHe+"
        , pBasisType       = "STO-3G"
        , pBasis           = [baseHe,baseH]
        , atomList         = [
                                  [0.0,0.0,1.4632]
                                 ,[0.000000,0.000000,0.000000]
                                ]
        }
    where
        sH = [(0.444635,0.168856),(0.535328,0.623913),(0.154329,3.42525)]
        sHe =[ (0.44463454,0.31364979),(0.53532814,1.15892300),(0.15432897,6.36242139)]
        baseH  = normGlobal`fmap` [CGF sH S]
        baseHe = normGlobal`fmap` [CGF sHe S]

        
project "CO" "STO-3G" = 
    PD {
          pLabel           = "CO"
        , pBasisType       = "STO-3G"
        , pBasis           = [baseO,baseC]
        , atomList         = [
                                  [0.0,0.0,0.0]
                                 ,[0.000000,0.000000,2.13539]
                                ]
        }
    where
        sO1 = [(0.444635,6.4436083),(0.535328,23.8088610),(0.154329,130.7093200)]
        sO2 = [(-0.09996723,5.0331513),(0.39951283,1.1695961),(0.7001154689,0.3803890)]
        pO =  [(0.15591627,5.0331513),(0.60768372,1.1695961),(0.39195739,0.3803890)]
        sC1 = [(0.1543289673,0.7161683735e2),(0.5353281423,0.1304509632e2),(0.4446345422,0.3530512160e1)]
        sC2 = [(-0.9996722919e-1,0.2941249355e1),(0.3995128261,0.6834830964),(0.7001154689,0.2222899159)]
        pC =  [(0.1559162750,0.2941249355e1),(0.6076837186,0.6834830964),(0.3919573931,0.2222899159)]
        baseC  = normGlobal`fmap` [CGF sC1 S, CGF sC2 S,CGF pC Px, CGF pC Py, CGF pC Pz ]
        baseO = normGlobal`fmap` [CGF sO1 S, CGF sO2 S, CGF pO Px, CGF pO Py, CGF pO Pz]
        
       
project "water" "STO-3G" =
    PD {
          pLabel           = "water"
        , pBasisType       = "STO-3G"
        , pBasis           = [baseH,baseO]
        , atomList         = [
                                  [0.000000,0.000000,0.000000]
                                , [0.000000,0.000000,1.814137]
                                , [1.710385,0.000000,-0.60471]
                                ]
        }
    where
        sO1 = [(0.444635,6.4436083),(0.535328,23.8088610),(0.154329,130.7093200)]
        sO2 = [(-0.09996723,5.0331513),(0.39951283,1.1695961),(0.7001154689,0.3803890)]
        pO = [(0.15591627,5.0331513),(0.60768372,1.1695961),(0.39195739,0.3803890)]
        sH = [(0.444635,0.168856),(0.535328,0.623913),(0.154329,3.42525)]
        baseH = normGlobal`fmap` [CGF sH S]
        baseO = normGlobal`fmap` [CGF sO1 S, CGF sO2 S, CGF pO Px, CGF pO Py, CGF pO Pz]

        
project "water" "6-31G*" =
    PD {
          pLabel           = "water"
        , pBasisType       = "6-31G*"
        , pBasis           = [baseH,baseO]
        , atomList         = [
                                  [0.000000,0.000000,0.000000]
                                , [0.000000,0.000000,1.81406]
                                , [1.710312,0.000000,-0.60469]
                                ]
        }
    where
        sO1 = [(0.1831074430e-2,0.5484671660e4),(0.1395017220e-1,0.8252349460e3),(0.6844507810e-1,0.1880469580e3),(0.2327143360,0.5296450000e2),(0.4701928980,0.1689757040e2),(0.3585208530,0.5799635340e1)]
        sO2 = [(-0.1107775495,0.1553961625e2),(-0.1480262627,0.3599933586e1),(0.1130767015e1,0.1013761750e1)]
        pO1 = [(0.7087426823e-1,0.1553961625e2),(0.3397528391,0.3599933586e1),(0.7271585773,0.1013761750e1)]
        sO3 = [(1.0,0.2700058226)]
        pO2 = [(1.0,0.2700058226)]
        dO  = [(1.0,0.8000000000)]
        sH1 = [(0.3349460434e-1,0.1873113696e2),(0.2347269535,0.2825394365e1),(0.8137573261,0.6401216923)]
        sH2 = [(1.0,0.1612777588)]
        baseH = normGlobal`fmap` [CGF sH1 S, CGF sH2 S]
        baseO = normGlobal`fmap` [CGF sO1 S, CGF sO2 S, CGF pO1 Px, CGF pO1 Py, CGF pO1 Pz, CGF sO3 S, CGF pO2 Px, CGF pO2 Py, CGF pO2 Pz, CGF dO Dxx, CGF dO Dyy, CGF dO Dzz, CGF dO Dxy,CGF dO Dxz, CGF dO Dyz]
