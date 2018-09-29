
# coding: utf-8

# # Helper methods
# 
# This notebook describes generic methods used throught the study.

# In[1]:

class PythonHelper(object):
    """ A generic Python helper class describing various helper methods.
    """

    # write to Jupyter console
    @staticmethod
    def writeToJupyterConsole(message):
        """  Writes text to Jupyter console.
             
             Args:
                 message: the message to write to the console.
        """
        import sys        
        nb_stdout = sys.stdout
        try:        
            sys.stdout = open('/dev/stdout', 'w')
            print(message)
        finally:
            sys.stdout = nb_stdout
    
    @staticmethod
    def getDecimal(text):
        """  Converts a given string to a Python Decimal.  It takes care to handle None values.
        
             Args:
                 text: the string to convert.
                 
             Returns:
                 A Decimal value representation of the input string or None.
        """
        from decimal import Decimal
        if text == None or text == 'None':
            return None
        else:
            return Decimal(text)


# In[ ]:



