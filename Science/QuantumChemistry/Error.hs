{-# Language DeriveDataTypeable #-}

{-|
Module: Science.QuantumChemistry.Error 
Description: Types for error reporting
Copyright: @2016 Felipe Zapata
-}


module Science.QuantumChemistry.Error where

-- ====================> Standard and third party libraries <========================
import Control.Exception 
import Data.Typeable


-- ===================> Data Types <===============================
-- | Error Reporting
data HSFOCKException = SCFError      -- ^ SCF related errors 
                     | ParseError -- ^ Input Output Related Errors
                     | KeyError      -- ^ Unkown key in Dictionary 
                       deriving (Show,Typeable)

instance Exception HSFOCKException

