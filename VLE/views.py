from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from django.forms import formset_factory
from .forms import FugacityVaporForm, ActivityForm
from .engine.fugacity import *
from .engine.chemsep_operation import get_chemical
from math import factorial


def vle_home(request):
    """VLE homepage"""

    return render(request, "vle_home.html", {"content":""})

def activity(request):
    """Bubble/dew point calcualtion page"""
    form = ActivityForm()
    return render(request, "activity.html", {"form":form})

def activity_result(request):

    #try:
    if request.method == 'GET':
        
        data = request.GET
        method = data['method']
        T = float(data['temperature'])
        P = float(data['pressure'])
        lang = int(data['compound_no'])

        chemicals = []; fractions = []; 
        
        for i in range(0,lang):
            chemicals.append(get_chemical(data['chemical'+str(i)]))
            fractions.append(float(data['fraction'+str(i)]))
        
        if len(chemicals)==1:
            fractions[0] = 1
        
        result = activities(chemicals, T, fractions, method)
        gammas = result[0]
        
        names=[]; cas=[]
        for i in range(0,len(chemicals)):
            names.append( chemicals[i].name )
            cas.append( chemicals[i].CAS )
        report = zip(names,cas,gammas)

        bips = {}
        if method != "Unifac" and method != "Dortmund":
            if method != "Ideal":
                for key in result[1].keys():
                    bips[ names[key[0]]+'-'+names[ key[1]] ] = result[1][key]
        
            else:
                bips = result[1]
        
        else:
            for key in result[1].keys():
                    bips["Groups "+str(key[0])+"-"+str(key[1]) ] = result[1][key]
            

        context = {"report":report,"bips":bips}
    #except:
        #context = {"report":"An error occured. Please check the parameters, methods and other inputs. For some property methods, there might be instable solution with entered parameters."}

    return render(request, "ajax/activity_result.html", context)

def PT_flash(request):
    """Pressure-temperature flash calcualtion page"""

    return render(request, "pt_flash.html", {"formset":""})


def fugacity(request):
    """Fugacity calcualtion page"""
    form = FugacityVaporForm()
    
    return render(request, "fugacity.html", {"form":form})


def fugacity_vapor_result(request):
    
    phis = ""
    try:
        if request.method == 'GET':
            
            data = request.GET
            method = data['method']
            kij_method = data['kij']
            T = float(data['temperature'])
            P = float(data['pressure'])
            lang = int(data['compound_no'])

            chemicals = []; fractions = []; 
            kij_input = None; kij_tune = None
            
            for i in range(0,lang):
                chemicals.append(get_chemical(data['chemical'+str(i)]))
                fractions.append(float(data['fraction'+str(i)]))
            
            if kij_method == "kij_input":
                kij_input = {}; kij_tune = None
                if lang == 2:
                    kij_input[(0,1)] = float(data['kij_input01'])
                
                elif lang == 3:
                    kij_input[(0,1)] = float(data['kij_input01'])
                    kij_input[(0,2)] = float(data['kij_input02'])
                    kij_input[(1,2)] = float(data['kij_input12'])
                
                elif lang == 4:
                    kij_input[(0,1)] = float(data['kij_input01'])
                    kij_input[(0,2)] = float(data['kij_input02'])
                    kij_input[(1,2)] = float(data['kij_input12'])
                    kij_input[(0,3)] = float(data['kij_input03'])
                    kij_input[(1,3)] = float(data['kij_input13'])
                    kij_input[(2,3)] = float(data['kij_input23'])

            elif kij_method == "tune":
                kij_input = None; kij_tune = {}
                
                if lang == 2:
                    kij_tune[(0,1)] = float(data['tune01'])
                
                elif lang == 3:
                    kij_tune[(0,1)] = float(data['tune01'])
                    kij_tune[(0,2)] = float(data['tune02'])
                    kij_tune[(1,2)] = float(data['tune12'])
                
                elif lang == 4:
                    kij_tune[(0,1)] = float(data['tune01'])
                    kij_tune[(0,2)] = float(data['tune02'])
                    kij_tune[(1,2)] = float(data['tune12'])
                    kij_tune[(0,3)] = float(data['tune03'])
                    kij_tune[(1,3)] = float(data['tune13'])
                    kij_tune[(2,3)] = float(data['tune23'])
                    
            if len(chemicals)==1:
                fractions[0] = 1
            
            result = fugacity_vapor(chemicals, T, P, fractions, method, kij_input, kij_tune)
            phis = result[0]
            
            names=[]; cas=[]
            for i in range(0,len(chemicals)):
                names.append( chemicals[i].name )
                cas.append( chemicals[i].CAS )
            report = zip(names,cas,phis)

            bips = {}
            if method != "Ideal":
                for key in result[1].keys():
                    bips[ names[key[0]]+'-'+names[ key[1]] ] = result[1][key]
            
            else:
                bips = result[1]
            
            context = {"report":report,"bips":bips}
    except:
        context = {"report":"An error occured. Please check the parameters, methods and other inputs. For some property methods, there might be instable solution with entered parameters."}

    return render(request, "ajax/fugacity_vapor_result.html", context)



