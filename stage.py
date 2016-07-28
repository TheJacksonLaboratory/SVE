#only have to instantiate this class and pass the wrapper
#to dynamically find and load the corresponding stage wrapper class
#can use pre registered stages from the db or read a JSON
import inspect
import stages

class Stage:
    #load up all wrappers in the import stages and use a poly reference to the wrapper
    #via the inspection module, this means that adding a new stage wrapper with its stage_id
    #will auto-magically load, bind and function correctly without using its name or constructor explicitly
    def __init__(self,wrapper='',dbc=None,retrieve=False,upload=True,params=None):
        #get all the wrapper classes
        ms,self.__link__ = inspect.getmembers(stages,inspect.ismodule),None   
        for m in ms: #[0] is the wrapper, [1] is the class definition
            if m[0]==wrapper:
                c = inspect.getmembers(m[1],inspect.isclass)[0][1]
                self.__link__ = c(wrapper,dbc,retrieve,upload,params)
        if self.__link__ is None:
            print('wrapper: %s not registered'%wrapper)
            self.params = {}
        else:
            print('using wrapper: %s'%wrapper)
            self.params = self.__link__.get_params()
            
    def get_params(self):
        return self.params
        
    def set_params(self,params) :
        self.__link__.set_params(params)
        
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return 0
    
    def run(self,run_id,inputs):
        return self.__link__.run(run_id,inputs)