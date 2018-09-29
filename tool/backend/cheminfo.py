class ChemInfo(object):

    @staticmethod
    def smilesToSVG(smiles):        
        """
            Args: 
                smiles : str in SMILES format.
                
            Returns:
                Str in SVG format.
        """
        
        import openbabel
        import re
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "svg")
        outMol = openbabel.OBMol()
        obConversion.ReadString(outMol, str(smiles))
        ans = obConversion.WriteString(outMol)
        
        ## Make the svg background transparent:
        ## replace fill="rgb(255,255,255)" with fill-opacity="0"
        ans = re.sub("fill=\"white\"", "fill-opacity=\"0\"", ans)
        
        return ans