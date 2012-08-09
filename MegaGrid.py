#from http://wxpython-users.1045709.n5.nabble.com/MegaGrid-Example-Shows-how-to-resize-virtual-grids-td2299061.html
import  wx, itertools, operator, math, numpy as np
import  wx.grid as  Grid

class MegaTable(Grid.PyGridTableBase):
    """
    A custom wx.Grid Table using user supplied data
    """
    def __init__(self, data, colnames, plugins):
        """data is a numpy array of the form:
        [rowname, col values]
        dictionary.get(colname, None) returns the data for column
        colname
        """
        # The base class must be initialized *first*
        Grid.PyGridTableBase.__init__(self)
        self.data = np.array(data)
        self.data = np.resize(len(data)/len(colnames),len(colnames)+1)
        self.hiddencols, self.colnames = [], []
#        self.colnames = colnames
        self.colIndexes = colnames
        for i in colnames:
            if i.find('.Hide') != -1:
                self.hiddencols.append(i)
            else:
                self.colnames.append(i)
        #self.hiddencols
        #self.colnames = colnames
        self.plugins = plugins or {}
        # XXX
        # we need to store the row length and column length to
        # see if the table has changed size
        self._rows = self.GetNumberRows()
        self._cols = self.GetNumberCols()

    def OnKey(self, event):
        print event.GetKeyCode()

    def GetNumberCols(self):
        return len(self.colIndexes)
    
    def GetVisibleCols(self):
        return len(self.colnames)

    def GetNumberRows(self):
        return len(self.data)

    def GetColLabelValue(self, col):
        return self.colnames[col]
    
    def GetHiddenColLabelValue(self, col):
        return self.hiddencols[col]
    
    def GetRowLabelValue(self, row):
        return "%s" % str(self.data[row,0])

    def GetValue(self, row, col):
        return str(self.data[row,col+1])
    
    def GetRawValue(self, row, col):
        return self.data[row,col+1]

    def SetValue(self, row, col, value):
        self.data[row,self.GetColLabelValue(col)] = value
        
    def ResetView(self, grid):
        """
        (Grid) -> Reset the grid view.   Call this to
        update the grid if rows and columns have been added or deleted
        """
        grid.BeginBatch()

        for current, new, delmsg, addmsg in [
            (self._rows, self.GetNumberRows(), Grid.GRIDTABLE_NOTIFY_ROWS_DELETED, Grid.GRIDTABLE_NOTIFY_ROWS_APPENDED),
            (self._cols, self.GetNumberCols(), Grid.GRIDTABLE_NOTIFY_COLS_DELETED, Grid.GRIDTABLE_NOTIFY_COLS_APPENDED),
        ]:
            if new < current:
                msg = Grid.GridTableMessage(self,delmsg,new,current-new)
                grid.ProcessTableMessage(msg)
            elif new > current:
                msg = Grid.GridTableMessage(self,addmsg,new-current)
                grid.ProcessTableMessage(msg)
                self.UpdateValues(grid)

        grid.EndBatch()

        self._rows = self.GetNumberRows()
        self._cols = self.GetNumberCols()
        # update the column rendering plugins
        self._updateColAttrs(grid)

        # update the scrollbars and the displayed part of the grid
        grid.AdjustScrollbars()
        grid.ForceRefresh()


    def UpdateValues(self, grid):
        """Update all displayed values"""
        # This sends an event to the grid table to update all of the values
        msg = Grid.GridTableMessage(self, Grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
        grid.ProcessTableMessage(msg)

    def _updateColAttrs(self, grid):
        """
        wx.Grid -> update the column attributes to add the
        appropriate renderer given the column name.  (renderers
        are stored in the self.plugins dictionary)

        Otherwise default to the default renderer.
        """
        col = 0

        for colname in self.colnames:
            attr = Grid.GridCellAttr()
            if colname in self.plugins:
                renderer = self.plugins[colname](self)

                if renderer.colSize:
                    grid.SetColSize(col, renderer.colSize)

                if renderer.rowSize:
                    grid.SetDefaultRowSize(renderer.rowSize)
                attr.SetReadOnly(True)
                attr.SetRenderer(renderer)

            grid.SetColAttr(col, attr)
            col += 1

    # ------------------------------------------------------
    # begin the added code to manipulate the table (non wx related)
    def AppendRow(self, inputData):
        inputData.insert(0,len(self.data))
        try:
            self.data = np.append(self.data,[inputData],axis=0)
        except ValueError:
            inputData[0]=0
            self.data = np.array([inputData])

    def DeleteCols(self, cols):
        """
        cols -> delete the columns from the dataset
        cols hold the column indices
        """
        # we'll cheat here and just remove the name from the
        # list of column names.  The data will remain but
        # it won't be shown
        self.data = np.delete(self.data,cols,axis=1)

        if not len(self.colnames):
            self.data = np.array([[]])

    def DeleteRows(self, rows):
        """
        rows -> delete the rows from the dataset
        rows hold the row indices
        """
        self.data = np.delete(self.data, rows, axis=0)
        
    def SetRowLabel(self, row, value):
        self.data[row,0] = value

    def SortColumn(self, col, reindex):
        """
        overridden
        col -> sort the data based on the column indexed by col
        """
        self.data = self.data[self.data[:,col].argsort()]
        if reindex: #changed so the row index is changed with sorting
            al = len(self.data)
            indices = np.linspace(0,al-1,al)
            self.data[:,0] = indices

    # end table manipulation code
    # ----------------------------------------------------------


# --------------------------------------------------------------------
# Sample wx.Grid renderers



class MegaFontRenderer(Grid.PyGridCellRenderer):
    """
    Changed from a generic font rendered to a specific one meant to display a grouped set of peptides
    """
    def __init__(self, table, color="blue", font="ARIAL", fontsize=8):
        """Render data in the specified color and font and fontsize"""
        Grid.PyGridCellRenderer.__init__(self)
        self.table = table
        self.color = color
        self.font = wx.Font(fontsize, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, font)
        self.selectedBrush = wx.Brush("blue", wx.SOLID)
        self.normalBrush = wx.Brush(wx.WHITE, wx.SOLID)
        self.colSize = None
        self.rowSize = 0

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        # Here we draw text in a grid cell using various fonts
        # and colors.  We have to set the clipping region on
        # the grid's DC, otherwise the text will spill over
        # to the next cell
        dc.SetClippingRect(rect)

        # clear the background
        dc.SetBackgroundMode(wx.SOLID)
        
        if isSelected:
            dc.SetBrush(wx.Brush(wx.BLUE, wx.SOLID))
            dc.SetPen(wx.Pen(wx.BLUE, 1, wx.SOLID))
        else:
            dc.SetBrush(wx.Brush(wx.WHITE, wx.SOLID))
            dc.SetPen(wx.Pen(wx.WHITE, 1, wx.SOLID))
        dc.DrawRectangleRect(rect)

        text = self.table.GetValue(row, col)
        dc.SetBackgroundMode(wx.SOLID)

        # change the text background based on whether the grid is selected
        # or not
        if isSelected:
            dc.SetBrush(self.selectedBrush)
            dc.SetTextBackground("blue")
        else:
            dc.SetBrush(self.normalBrush)
            dc.SetTextBackground("white")

        dc.SetTextForeground(self.color)
        dc.SetFont(self.font)
        dc.DrawText(text, rect.x+1, rect.y+1)

        # Okay, now for the advanced class :)
        # Let's add three dots "..."
        # to indicate that that there is more text to be read
        # when the text is larger than the grid cell

        width, height = dc.GetTextExtent(text)
        
        if width > rect.width-2:
            width, height = dc.GetTextExtent("...")
            x = rect.x+1 + rect.width-2 - width
            dc.DrawRectangle(x, rect.y+1, width+1, height)
            dc.DrawText("...", x, rect.y+1)

        dc.DestroyClippingRegion()


# --------------------------------------------------------------------
# Sample Grid using a specialized table and renderers that can
# be plugged in based on column names

class MegaGrid(Grid.Grid):
    def __init__(self, parent, data, colnames, plugins=None):
        """parent, data, colnames, plugins=None
        Initialize a grid using the data defined in data and colnames
        (see MegaTable for a description of the data format)
        plugins is a dictionary of columnName -> column renderers.
        """

        # The base class must be initialized *first*
        Grid.Grid.__init__(self, parent, -1)
        self._table = MegaTable(data, colnames, plugins)
        self.SetTable(self._table)
        self._plugins = plugins
        self.cellSelection = []
        self.parent = parent

        self.Bind(Grid.EVT_GRID_LABEL_RIGHT_CLICK, self.OnLabelRightClicked)
        self.Bind(Grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
        self.Bind(Grid.EVT_GRID_SELECT_CELL, self.OnCellSelect)
        self.Bind(Grid.EVT_GRID_RANGE_SELECT, self.OnRangeSelect)
        self.Bind(wx.EVT_KEY_UP, self.OnKey)

    def Reset(self):
        """reset the view based on the data in the table.  Call
        this when rows are added or destroyed"""
        self._table.ResetView(self)
        
    def OnKey(self, event):
        key = event.GetKeyCode()
        if key == 67:
            if event.ControlDown():
                #ctrl+c pressed
                if (len(self.cellSelection) == 1):
                    entry = self.GetCellValue(self.cellSelection[0][0], self.cellSelection[0][1])
                else:
                    #sort by row
                    self.cellSelection.sort(key=lambda x:x[0])
                    lastRow = self.cellSelection[0][0]
                    entry = ""
                    for row, j in itertools.groupby(self.cellSelection, key=lambda x:x[0]):
                        cols = [k[1] for k in j]
                        cols.sort()
                        for col in cols:
                            if row != lastRow:
                                lastRow = row
                                entry = entry[:-1] #remove trailing tab
                                entry += "\n"+self.GetCellValue(row,col)+"\t"
                            else:
                                entry += self.GetCellValue(row,col)+"\t"
                self.toClipboard(entry)
        event.Skip()

    

    def OnRangeSelect(self, event):
        tlIndex = (event.GetTopRow(), event.GetLeftCol())
        brIndex = (event.GetBottomRow(), event.GetRightCol())
        if event.Selecting():
            #adding to the list...
            for rIndex in range( tlIndex[0], brIndex[0]+1):
                for cIndex in range(tlIndex[1], brIndex[1]+1):
                    selCell = tuple((rIndex,cIndex))
                    if selCell not in self.cellSelection:
                        self.cellSelection.append( selCell )
        else:
            if event.ControlDown():
                # removal from list
                for rIndex in range( tlIndex[0], tlIndex[1]+1):
                    for cIndex in range(brIndex[0], brIndex[1]+1):
                        selCell = tuple((rIndex,cIndex))
                        if selCell in self.cellSelection:
                            self.cellSelection.remove( selCell )
            else:
                pass
        event.Skip()
            
    
    def OnCellSelect(self, event):
        self.cellSelection = [(event.GetRow(), event.GetCol())]
        event.Skip()
    
    def SetValue(self, row, col, entry):
        self._table.SetValue(row,col,entry)
        
    def GetRowValue(self, row):
        return self._table.GetRowLabelValue(row)
    
    def RemoveRows(self, rows):
        self._table.DeleteRows(rows)
        
    def AppendRow(self, data, reset):
       self._table.AppendRow(data)
       if reset:
           self.Reset()
       
    def InsertRow(self, rowindex, data):
        self._table.InsertRow(rowindex, data)
       
    def OnCellRightClick(self, event):
       row, col = event.GetRow(), event.GetCol()
       self.cellPopup(row, col, event)
       
    def toClipboard(self, text):
       clipdata = wx.TextDataObject()
       clipdata.SetText(text)
       wx.TheClipboard.Open()
       wx.TheClipboard.SetData(clipdata)
       wx.TheClipboard.Close()
    
    def OnLabelRightClicked(self, evt):
        # Did we click on a row or a column?
        row, col = evt.GetRow(), evt.GetCol()
        if row == -1: self.colPopup(col, evt)

    def colPopup(self, col, evt):
        """(col, evt) -> display a popup menu when a column label is
        right clicked"""
        x = self.GetColSize(col)/2
        menu = wx.Menu()
        id1 = wx.NewId()
        sortID = wx.NewId()

        xo, yo = evt.GetPosition()
        self.SelectCol(col)
        cols = self.GetSelectedCols()
        self.Refresh()
#        menu.Append(id1, "Delete Col(s)")
        menu.Append(sortID, "Sort Column")

        def sort(event, self=self, col=col):
            self._table.SortColumn(col, False)
            self.Reset()

        if len(cols) == 1:
            self.Bind(wx.EVT_MENU, sort, id=sortID)

        self.PopupMenu(menu)
        menu.Destroy()
        return

        
class MegaFontRendererFactory:
    def __init__(self, color, font, fontsize):
        """
        (color, font, fontsize) -> set of a factory to generate
        renderers when called.
        func = MegaFontRenderFactory(color, font, fontsize)
        renderer = func(table)
        """
        self.color = color
        self.font = font
        self.fontsize = fontsize

    def __call__(self, table):
        return MegaFontRenderer(table, self.color, self.font, self.fontsize)
#
#
##---------------------------------------------------------------------------
#
#class TestFrame(wx.Frame):
#    def __init__(self, parent, plugins={"This":MegaFontRendererFactory("red", "ARIAL", 8),
#                                        "Test":MegaFontRendererFactory("orange", "TIMES", 24),}):
#        wx.Frame.__init__(self, parent, -1,
#                         "Test Frame", size=(640,480))
#
#        grid = MegaGrid(self, data, colnames, plugins)
#        grid.Reset()
#

