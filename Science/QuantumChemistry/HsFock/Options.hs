
-- The thing...
-- A system shell for the HaskellFock SCF Project 
-- @2012,2013 Angel Alvarez Adhesive tape
-- @2012,2013 Felipe Zapata core SCF machinery 

module Science.QuantumChemistry.HsFock.Options where

import Text.Show.Functions



-- Record for storing cmdline options
data Options = Options
 { 
--  optDump        :: Bool
   optModules     :: [(String, Options-> IO() ) ]
 , optMode        :: Maybe (Options->IO())
 , optVerbose     :: Bool
 , optShowVersion :: Bool
 , optOutput      :: Maybe FilePath
 , optDataDir     :: Maybe FilePath
 , optInput       :: [FilePath]
 } deriving (Show)
