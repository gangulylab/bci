import pygame as pg
import numpy as np
import time
import os
import csv

class UI(object):
	def __init__(self, log, updateRate):
		self.screen = pg.display.get_surface()
		self.done = False 
		self.mode = 0
		self.makeButtons()
		self.newTime = time.time()
		self.modeList = ['Translation', 'Orientation', 'Gripper']
		self.colorList = [pg.Color('dodgerblue'), pg.Color('red'), pg.Color('green')]
		self.font = pg.font.Font(None, 40)
		self.updateText()
		self.buttonState = -1
		self.state = -1
		self.prevbuttonState = -1
		self.log = log
		self.updateRate = updateRate
		

		if self.log:
			dataDir = time.strftime("%Y%m%d")

			if not os.path.exists(dataDir):
				os.makedirs(dataDir)
			trialInd = len(os.listdir(dataDir))
			self.fn = dataDir + "/data00" + str(trialInd) + ".csv"
			self.startLogger(self.fn)

	def makeButtons(self):
		self.xl = []
		self.xu = []
		self.yl = []
		self.yu = []
		button = []
		rect = []
		self.screen.fill((0,0,0))

		for i in range(9):
			button.append(pg.image.load('img/c' + str(self.mode) + '_' + str(i+1) + '.png').convert())
			rect.append(button[i].get_rect())
			row = i // 3
			col = np.mod(i,3)
			rect[i].center = 100 + 150*col, 400 - row*150
			self.xl.append(rect[i].x)
			self.xu.append(rect[i].x + rect[i].w)
			self.yl.append(rect[i].y)
			self.yu.append(rect[i].y + rect[i].h)
			self.screen.blit(button[i], rect[i])
			self.x = 0
			self.y = 0
			pg.draw.rect(self.screen, (0,0,0), rect[i], 1)

		pg.display.update()

	def updateText(self):

		self.makeButtons

		pgTxt1 = 'Mode ' + str(self.mode + 1) + ': ' + self.modeList[self.mode]
		color = self.colorList[self.mode]

		self.screen.fill((0,0,0), (0,0,300, 50))  
		self.txt_surface1 = self.font.render(pgTxt1, True, color)
		self.screen.blit(self.txt_surface1, (10, 10))
		pg.display.flip()

	def update(self):	

		for event in pg.event.get():
			if event.type == pg.QUIT:
				self.done = True

		m = pg.mouse.get_pos()
		self.x = m[0]
		self.y = m[1]

		if time.time() - self.newTime > self.updateRate :
			self.newTime = time.time()

			self.prevbuttonState = self.buttonState
			self.buttonState = -1
			self.state = -1
			
			for i in range(9):
				if  (self.xl[i]  <= m[0] <=  self.xu[i]) & (self.yl[i]  <= m[1] <=  self.yu[i]):
					self.buttonState = i+1

			if self.buttonState == 2:
					self.state = 4
			elif self.buttonState == 4:
				self.state = 3
			elif self.buttonState == 6:
				self.state = 1
			elif self.buttonState == 8:
				self.state = 2

			# if self.buttonState == 3 and self.prevbuttonState != 3:
			# 	if self.mode < 1:
			# 		self.mode = self.mode + 1
			# 	else:
			# 		self.mode = 0
			# 	self.updateText()	
		if self.log:
			self.updateLogger()

	def startLogger(self, fn):
			self.file  = open(fn, "w", newline = '')
			self.fileObj = csv.writer(self.file)
			self.logOpen = True
			self.logT = time.time()

	def updateLogger(self):

		if time.time() - self.logT > 0.01:
			self.logT = time.time()
			line = [time.time(), self.x, self.y, self.buttonState, self.mode]
			self.fileObj.writerow(line)

	def closeLogger(self):
		self.file.close()
		self.logOpen = False	
			

if __name__ == "__main__":
	pg.init()
	pg.display.set_mode((500,500))
	pg.display.set_caption("Control Interface")
	runUI = UI(1, 0.3)

	while not runUI.done:
		for event in pg.event.get():
			if event.type == pg.QUIT:
				runUI.closeLogger()
				runUI.done = True
		runUI.update()
	pg.quit()