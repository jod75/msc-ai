
# coding: utf-8

# # HDFS Helper class
# 
# This helper requires `hdfs` module.  To install use:
# 
# ```bash
# sudo pip install hdfs
# ```

# In[ ]:

class HDFSHelper(object):       
    
    @staticmethod
    def fileExists(server, filename):   
        """ Checks is a file exists in hdfs.
        
            Args:
                server: hdfs server and port, example: "http://hadoop1:50070".
                filename: path of file to check.
                
            Returns:
                True: if file exists, False otherwise.
        """
        
        from hdfs import InsecureClient
        
        client = InsecureClient(server)
        return client.content(filename, strict=False) != None
    
    @staticmethod
    def putFile(server, source, destination):
        """ Uploads a file to HDFS.
        
            Args:
                server: hdfs server and port, example: "http://hadoop1:50070".
                source: localpath of file to upload. 
                destination: remote file, inclusing path.
        """
        from hdfs import InsecureClient
        
        client = InsecureClient(server)
        return client.upload(hdfs_path=destination, local_path=source)
        


# In[ ]:



