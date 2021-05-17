#!/usr/bin/env python
# coding=UTF-8
'''
Author: Howie Hong(希理)
Date: 2020-01-10 09:41:39
LastEditors: Howie Hong(希理)
LastEditTime: 2021-01-15 08:27:24
Description: 
'''
from matplotlib.pyplot import xlim,ylim
import numpy,multiprocessing
from scipy.constants import G
import matplotlib.pyplot as plt
from os import _exit as exit

def magnitude(vector):
    '''
    description: calculate the magnitude of a numpy vector

    parameter: vector{numpy.array}

    return: magnitude of the numpy vector{float}
    '''
    return numpy.sqrt(vector.dot(vector))

class celestial_body:
    '''
    celestial_body class
    '''
    def __init__(self,name,mass,location,velocity,color):
        '''
        parameter:
            name{str}: name of the celestial body
            mass{float}: mass of the celestial body, in kg
            location{2D float numpy.array}: relative location of the celestial body, in m
            velocity{numpy.array} :relative velocity of the celestial body, in m/s
        '''
        self.name = name
        self.mass = mass
        self.location_raw = multiprocessing.Array('d',[location[0],location[1]])
        
        self.velocity_raw = multiprocessing.Array('d',[velocity[0],velocity[1]])
        
        self.color = color
    
    def __sync(self,synchronizer):
        '''
        description:
            stop and wait for other processes

        parameter: 
            synchronizer{sync Barrier} 

        return: 
        '''
        synchronizer.wait()
        if synchronizer.broken == True:
            synchronizer.reset()

    def orbit(self,other_celestial_body_list,time,synchronizer,debug=False):
        location = numpy.array([self.location_raw[0],self.location_raw[1]])
        velocity = numpy.array([self.velocity_raw[0],self.velocity_raw[1]])
        while True:
            '''
            description: 

            parameter:
                other_celestial_body_list: a list of celestial_body
                time{float}: bigger than 0, the smaller it is, the more precise the result is
                synchronizer{sync Barrier}: to sync between different processes
                debug{bool}: enable and disable debug, must be enable with gui

            return: 
            '''    
            c = 0.0 #a variable dependent on radius from other bodies and mass of other bodies
            for other_celestial_body in other_celestial_body_list:
                #convert synced array to numpy.array, and get the radius
                #r = numpy.array([other_celestial_body.location_raw[0],other_celestial_body.location_raw[1]]) - location
                r = other_celestial_body.location_raw - location
                c += (other_celestial_body.mass / numpy.dot(r,r)) * r/magnitude(r)

            location += velocity*time + G*c*(time**2) #G is imported from scipy
            velocity += G*c*time           
            
            self.__sync(synchronizer)
            self.location_raw[0] = location[0]
            self.location_raw[1] = location[1]
            
            if debug:
                self.velocity_raw[0] = velocity[0]
                self.velocity_raw[1] = velocity[1]
            
            self.__sync(synchronizer)


class process_manager:
    '''
    create and start processes
    '''      
    def __init__(self,celestial_body_list,time,trail=True,tail_length=100,debug=False):
        '''
        description: 

        parameter: 

        return: 
        '''
        self.celestial_body_list = celestial_body_list
        self.synchronizer = multiprocessing.Barrier(len(celestial_body_list))
        self.debug = debug
        self.trail = trail
        self.tail_length = tail_length
        self.time = time
        self.processes = []
        for n in range(len(self.celestial_body_list)):
            self.processes.append(multiprocessing.Process(target = celestial_body_list[n].orbit, \
                args= (celestial_body_list[:n]+celestial_body_list[n+1:],time,self.synchronizer,self.debug) ))

    def __resize(self,x_min,x_max,y_min,y_max):
        self.center = ((x_max+x_min)/2,(y_max+y_min)/2)
        self.size = max((x_max-x_min)/2,(y_max-y_min)/2)*1.1

    def __set_axis(self):
        xlim(self.center[0]-self.size,self.center[0]+self.size)
        ylim(self.center[1]-self.size,self.center[1]+self.size)
        plt.gca().set_aspect(1)

    def __gui(self):
        plt.ion()
        fig = plt.figure()
        fig.canvas.mpl_connect('close_event', self.__handle_close) #close the program when the window is close
        celestial_body_point_list = []
        x_min = x_max = y_min = y_max = 0
        for n in range(len(self.celestial_body_list)):
            celestial_body_point_list.append([[],[]])
        while True:
            plt.clf()
            for n in range(len(self.celestial_body_list)):
                if self.trail:
                    while len(celestial_body_point_list[n][0]) > self.tail_length:
                        celestial_body_point_list[n][0].pop(0)
                        celestial_body_point_list[n][1].pop(0)
                    celestial_body_point_list[n][0].append(self.celestial_body_list[n].location_raw[0])
                    celestial_body_point_list[n][1].append(self.celestial_body_list[n].location_raw[1])

                    plt.plot(celestial_body_point_list[n][0],celestial_body_point_list[n][1],'-',color = self.celestial_body_list[n].color)
                
                if self.celestial_body_list[n].location_raw[0] < x_min:
                    x_min = self.celestial_body_list[n].location_raw[0]
                    self.__resize(x_min,x_max,y_min,y_max)
                elif self.celestial_body_list[n].location_raw[0] > x_max:
                    x_max = self.celestial_body_list[n].location_raw[0]
                    self.__resize(x_min,x_max,y_min,y_max)
                if self.celestial_body_list[n].location_raw[1] < y_min:
                    y_min = self.celestial_body_list[n].location_raw[1]
                    self.__resize(x_min,x_max,y_min,y_max)
                elif self.celestial_body_list[n].location_raw[1] > y_max:
                    y_max = self.celestial_body_list[n].location_raw[1]
                    self.__resize(x_min,x_max,y_min,y_max)
                    
                plt.plot(self.celestial_body_list[n].location_raw[0],self.celestial_body_list[n].location_raw[1],'.',color = self.celestial_body_list[n].color)
            self.__set_axis()
            plt.pause(0.001)

            #the momentum is just a referencve value, it is not thread safety
            if self.debug:
                momentum = 0
                for celestial_body in self.celestial_body_list:
                    print(celestial_body.name,'location:',celestial_body.location_raw[0],celestial_body.location_raw[1])
                    momentum += celestial_body.mass * numpy.array([celestial_body.velocity_raw[0],celestial_body.velocity_raw[1]])
                print('total momentum:',momentum)
                print('time:',self.time)
                print()
            
    def start(self):
        for process in self.processes:
            process.start()
        self.__gui()
    
    def __handle_close(self,event):
        '''
        https://matplotlib.org/gallery/event_handling/close_event.html
        '''
        for process in self.processes:
            process.terminate()
        exit(0)

def main():
    celestial_body_list = []

    #PHYS 206 q1
    #celestial_body_list.append(celestial_body('Earth',5.965*10**24,(0.0,0.0),(0.0,0.0),'b'))
    #celestial_body_list.append(celestial_body('satellite',1,(6.4*10**6,0.0) ,(0.0,7890),'g'))#low-Earth orbit
    #celestial_body_list.append(celestial_body('satellite',1,(6.4*10**6,0.0) ,(0.0,10408,8),'g'))#elliptical orbit
    #celestial_body_list.append(celestial_body('satellite',1,(-4.22*10**7,0.0) ,(0.0,-3072.8),'g'))#geosynchronous orbit
    
    #PHYS 206 q2
    #celestial_body_list.append(celestial_body('celestial_body5',10**20,(0.0,0.0) ,(0.0,0.0),'g'))
    #celestial_body_list.append(celestial_body('celestial_body6',10**1,(10**5,0.0) ,(0.0,-1000),'b'))#The second last arg is velocity


    #test group 1
    #celestial_body_list.append(celestial_body('Sun',1.98191*10**30,(0.0,0.0),(0.0,0.0),'r'))
    #celestial_body_list.append(celestial_body('Earth',5.965*10**24,(1.521*10**11,0.0),(0.0,-29300.0),'b'))

    #test group 2
    celestial_body_list.append(celestial_body('celestial_body1',10**11,(0.0,100.0),(0.1,0.0),'r'))
    celestial_body_list.append(celestial_body('celestial_body2',10**11,(0.0,-100.0) ,(-0.1,0.0),'b'))
    celestial_body_list.append(celestial_body('celestial_body3',10**10,(100.0,0.0) ,(0.0,-0.3),'y'))
    celestial_body_list.append(celestial_body('celestial_body4',10**10,(-100.0,0.0) ,(0.0,0.3),'c'))
    celestial_body_list.append(celestial_body('celestial_body5',10**10,(0.0,0.0) ,(0.0,0.0),'g'))
    celestial_body_list.append(celestial_body('celestial_body6',10**10,(0.0,50.0) ,(0.15,0.0),'g'))
    celestial_body_list.append(celestial_body('celestial_body7',10**10,(0.0,-50.0) ,(-0.15,0.0),'g'))
    celestial_body_list.append(celestial_body('celestial_body8',10**10,(50.0,0.0) ,(0.0,-0.15),'g'))
    celestial_body_list.append(celestial_body('celestial_body9',10**10,(-50.0,0.0) ,(0.0,0.15),'g'))

    #celestial_body_list.append(celestial_body('celestial_body10',10**10,(100.0,100.0) ,(0.25,-0.25),'g'))
    #celestial_body_list.append(celestial_body('celestial_body11',10**10,(-100.0,100.0) ,(0.25,0.25),'g'))
    #celestial_body_list.append(celestial_body('celestial_body12',10**10,(100.0,-100.0) ,(-0.25,-0.25),'g'))
    #celestial_body_list.append(celestial_body('celestial_body13',10**10,(-100.0,-100.0) ,(-0.25,0.25),'g'))

    #celestial_body_list.append(celestial_body('celestial_body15',10**10,(50.0,50.0) ,(0.2,-0.2),'g'))
    #celestial_body_list.append(celestial_body('celestial_body16',10**10,(-50.0,50.0) ,(0.2,0.2),'g'))
    #celestial_body_list.append(celestial_body('celestial_body17',10**10,(50.0,-50.0) ,(-0.2,-0.2),'g'))
    #celestial_body_list.append(celestial_body('celestial_body18',10**10,(-50.0,-50.0) ,(-0.2,0.2),'g'))

    #test group 4
    #celestial_body_list.append(celestial_body('celestial_body1',10**11,(0.0,100.0),(0.15,0.0),'r'))
    #celestial_body_list.append(celestial_body('celestial_body3',10**10,(100.0,0.0) ,(0.0,-0.15),'y'))

    gui = process_manager(celestial_body_list,time=0.01,trail=True,tail_length=10**4,debug=False)
    gui.start()

if __name__ == "__main__":
    main()

