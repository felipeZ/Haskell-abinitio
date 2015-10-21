{-# LANGUAGE OverloadedStrings #-}

module Science.QuantumChemistry.BasisSet.FetchBasis where

-- ==================> Standard and third party libraries <=================
import           Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B
import           Network.Curl.Download (openURI)

-- =================> Internal Modules <======================

-- | Fetch the basis set from https://bse.pnl.gov/bse/portal
-- fetchBasis :: 

