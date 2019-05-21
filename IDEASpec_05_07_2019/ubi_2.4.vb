Option Strict Off
Option Explicit On
Imports VB = Microsoft.VisualBasic
Imports System.Windows.Forms
Imports System.Drawing
Imports Microsoft.VisualBasic.Compatibility
Imports ZedGraph
Imports System.IO

'vnew

Friend Class DMK_v1
    Inherits System.Windows.Forms.Form ' must be first (stupid VB...)
    Dim invert_raw As Integer = 1 'placve -1 here to invert the raw data (i.e. raw data will be multipolied by this)
    Dim version As String = "WASP version"
    Dim plot_data_point As Integer = True

    Dim program_directory As String = Environment.CurrentDirectory
    'Dim library_directory As String = program_directory & "\lib\"
    Dim rootdata_directory As String = program_directory & "\rootdata\"

    Dim Gain_Volts(2, 4, 10) As Single

    Dim Baseline_Y(5000) As Single
    Dim Auto_File_Name As String
    Dim auto_run As Integer = False
    Const Blue_Filter As Short = 110
    Const IR_Filter As Short = 165
    Dim Halt_Script As Short
    Dim Script() As String
    Dim record_files As Boolean = False
    Dim Script_Counter As Short
    Dim C_Script As String
    Public Graph_Name() As String = {"0", "1", "2", "3", "4", "5", "6"}
    Dim fluorescence_shutter_status As Integer
    Dim Shutter_delay As Integer = 1.6
    Const BoardNum As Short = 0 ' Board number
    Dim Trace_protocol() As Integer = {0, 0, 0} ' associates a protocol to a trace, such that the trace will automatically pick it's protocol when m_trace is called (deprecated)
    Dim Xe_intensity_value As Single
    Dim M_Primary_Gain(,) As Short = {{0, 0, 0}, {0, 0, 0}}
    Dim M_Reference_Gain(,) As Short = {{0, 0, 0}, {0, 0, 0}}

    Dim m_set_auto_gain(,) As Short = {{0, 0, 0}, {0, 0, 0}}
    Dim m_take_data(,) As Integer

    Dim Pre_Delay(,) As String = {{"0", "0", "0"}, {"0", "0", "0"}}

    Dim script_file As String
    Dim s_Plot_File_Name(10) As String
    Dim ref_mode As Integer
    Dim s_Plot_File_Type(10) As Integer
    Dim D_Out_Mask As Integer
    Dim laser As Integer = 0
    Dim Current_Wheel_Position As Single
    Dim view_laser As Integer
    Dim servo_position(8, 30) As Integer
    Dim current_servo, number_servo_positions, current_servo_position As Integer
    Dim sub_name As String
    Dim Sub_Return_Pointer(36) As Integer 'where to return after a subroutine
    Dim Sub_Return_Level As Integer
    Dim number_variables As Integer
    Dim user_variable(100, 100) As String
    Dim user_variable_name(1000) As String
    Dim user_variable_index(100) As Integer
    Dim Base_File_Name, Base_Base_File_Name As String
    Dim run_another_script As Integer = False
    Dim dice_length As Integer
    Dim take_notes As Integer
    Dim stir_status As Integer = 0
    Dim Current_Trace As Short
    Dim pre_pulse_light As Short
    Dim Max_Number_Lights As Integer = 32
    'from pblast"
    Dim stemp As Single
    'Dim Count_time_num$
    Dim temp1 As Short
    Dim temp2 As Short
    Dim temp3 As Short
    Dim PB_multiplier As Single
    Dim PB_freq As Double
    Dim PB_Numm As Object
    Dim PB_pulselen As Single
    Dim ptb1, ptb0, ptb2 As Short 'changed from object
    Dim ptb3 As Short
    Dim PB_cycles As Double
    Dim pb_flags As Integer
    Dim tempi As Integer
    Dim ccc(10) As Integer
    Dim tttemp As Short
    Dim PBStatt As Integer
    Dim pb_port As Short
    Dim max_measuring_pulse_duration As Single = 0.0002
    Dim nmp As Integer
    Dim new_firgelli As Integer = False

    Dim gain_slop As Single = 5.0

    Dim number_of_channels_to_scan As Integer

    Dim cblaster As New ChrisBlasterIO("Nexys2") ' object for communicating and using the ChrisBlaster FPGA board
    Dim protocol_menu() As Protocol
    ' Saturation Pulse button protocol
    Dim saturation_test As Protocol
    'Dim geoffs_detector_position As Integer = -9999
    'Dim absorbance_target_position As Integer = 4300
    'Dim fluorescence_target_position As Integer = 2900
    'Dim df_target_position As Integer = 1500

    'Dim WithEvents sp As New SerialPort
    Dim spec_mode As Integer = False 'False
    Dim start_lambda As Single = 525
    Dim delta_lambda As Single = 1
    Dim mono_BoardNum As Integer = 1 ' 1608fs board for monochromator control
    Dim mono_D_Out_Mask As Integer = 0 'mask for the digital IO for the monochromator board
    Dim steps_per_nanometer As Single = 96 / 5
    Dim current_lambda As Single = 0 'holds the value of the current monochromator position
    Dim x_mode As Integer = 0

    Dim list_of_base_additions As New System.Collections.Generic.List(Of String)



    Public Sub dmkversion1_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MyBase.Load
        'This sub places a default image in the pictureboxes to populate the pixels

        Dim XSize As Integer = 320
        Dim YSize As Integer = 240
        'Dim MyPic As String = "C:\Documents and Settings\All Users\Documents\My Pictures\Sample Pictures\hangar1.jpg"
        'Dim MyImage As New Bitmap(MyPic)
        'Dim s_Plot_File_name(10) As String

        ' Stretches the image to fit the pictureBoxes. 
        Dim i As Integer

        For i = 1 To 6
            Graph_Name(i) = "GRAPH: " & i
        Next i
        movelin.Visible = True

        '        Button8.Visible = False
        '        Button9.Visible = False
        '        TrackBar1.Visible = False
        '        Label6.Visible = False

        view_laser = False

        'set up variable to hold autogain stuff
        'iv = VB.Len(Script(Script_Counter)) - 1
        '       number_variables = 1
        '       user_variable_name(number_variables) = "auto_gain"
        '       user_variable_index(number_variables) = 1
        '        Jabber_jabber.List1.Items.Add(user_variable_name(number_variables))
        
        'this is for the monochromator on Charley
        'Dim BaudRates() As String = {"300", "1200", "2400", "4800", "9600", "14400", "19200", "28800", "38400", "57600", "115200"}
        'cmbBaud.Items.AddRange(BaudRates)
        'cmbBaud.SelectedIndex = 4
        'Try
        ' GetSerialPortNames()
        ' cmbPort.SelectedIndex = 0
        ' Catch
        ' MsgBox("No ports connected.")
        ' End Try
        Jabber_jabber.Show()
        'jabber_box.
        Jabber_jabber.List1.Items.Add("Hello!")
    End Sub

    Private Sub parsenum(ByRef num As String)
        'back up copy of code inserted into Private Sub P_blast

        Freq = 10 ^ 8
        pulselen = 1 / Freq

        StdDelay = 5 'standard 5 cycle delay with this board (or 50 nsec)

        num = Trim(num)

        If VB.Right(num, 1) = "u" Then
            multiplier = 0.000001
        ElseIf VB.Right(num, 1) = "m" Then
            multiplier = 0.001
        ElseIf VB.Right(num, 1) = "n" Then
            multiplier = 0.000000001
        Else
            multiplier = 1
        End If

        'seconds = multiplier * numm

        numm = Val(num)
        Time_Interval_in_Seconds = numm * multiplier

        'DMK_v1.List1.AddItem "n " & num$ & " TIS " & Time_Interval_in_Seconds & " m " & multiplier

        System.Windows.Forms.Application.DoEvents()

        cycles = Int(multiplier * numm / pulselen) - StdDelay

        'DMK_v1.List1.AddItem "cycles " & cycles

        System.Windows.Forms.Application.DoEvents()

        ptb0 = cycles And &HFFS
        cycles = (cycles - ptb0) / 256

        ptb1 = cycles And &HFFS
        cycles = (cycles - ptb1) / 256

        ptb2 = cycles And &HFFS
        cycles = (cycles - ptb2) / 256

        ptb3 = cycles And &HFFS
        'DMK_v1.List1.AddItem "0: " & ptb0 & " 1: " & ptb1 & " 2: " & ptb2 & " 3: " & ptb3

    End Sub

    Private Sub SaveToChrisBlaster(ByRef Current_Protocol As Short, ByVal memory_index As Integer, ByRef Number_Loops() As Short, _
                                   ByRef Intensity(,) As Short, ByRef Number_Pulses(,) As Short, ByRef Number_Measuring_Lights() As Short, _
                                   ByRef Measuring_Light(,) As Short, ByRef Choke_actinic(,) As Short, ByRef Measuring_Interval(,,) As String, _
                                   ByRef Primary_Gain(,) As Short, ByRef reference_gain(,) As Short, ByRef Measuring_Pulse_Duration(,) As String, _
                                   ByRef Gain() As Integer, ByRef In_Channel() As Short, _
                                   ByRef Ref_channel_number() As Short, _
                                   ByRef points As Integer, ByRef data_time() As Single, _
                                   ByRef Xe_Flash(,) As Short, ByRef q_switch(,) As Short, ByVal M_Far_Red(,) As Short, _
                                   ByVal M_Blue_actinic(,) As Short, ByVal pre_pulse(,) As Short, ByVal pre_pulse_time(,) As String, _
                                   ByVal Pre_Delay(,) As String, ByVal m_take_data(,) As Short)

        'Dim ulstat As Short


        Dim PP_Line As Short = 0
        'Dim i, j As Object
        Dim i As Object
        'Dim ii As Integer
        'Dim k As Short
        'Dim AD_Trigger As Short
        Dim Loop_Return_Index() As Short
        'Dim Options As Integer
        'Dim ss As Short
        'Dim CurCount, CurIndex As Integer
        Dim Time_Accumulator As Single = 0
        Dim Measuring_Pulse_Number As Integer = 0
        'Dim Timed_Out As Short
        Dim Time_Out_Time As Single = 1
        Dim Each_Measuring_Pulse_Interval() As Single
        Dim Current_Measuring_Point As Integer = 0
        'Dim Give_Saturating_Pulse As Short
        'Dim SH_State As Integer

        ReDim Loop_Return_Index(Number_Loops(Current_Protocol))
        ReDim Each_Measuring_Pulse_Interval(Number_Measuring_Lights(Current_Protocol))
        'Dim pbstat As Integer


        points = 0

        For i = 1 To Number_Loops(Current_Protocol)
            If m_take_data(Current_Protocol, i) > 0 Then
                points = points + (Number_Measuring_Lights(Current_Protocol) * Number_Pulses(Current_Protocol, i))
            End If
        Next i





        Dim expected_run_time As Double
        ' Use the ChrisBlaster to run the protocol
        Dim program As Protocol = New Protocol()

        ' The protocol object holds all of the relevent values of a protocol in one easy package for use by the ChrisBlasterIO object
        ' to handle the nuts & bolts of programming the hardware
        program.Adc_Gain = Gain(Current_Protocol)

        'program.Baseline_End = 
        'program.Baseline_Start = 
        program.In_Channel = In_Channel(Current_Protocol)
        program.laser = laser
        For n_light As Integer = 1 To Number_Measuring_Lights(Current_Protocol)
            'program.Measuring_Pulse_Duration(n_light) = Measuring_Pulse_Duration(Current_Protocol, n_light)

        Next
        program.Number_Loops = Number_Loops(Current_Protocol)
        program.Number_Measuring_Lights = Number_Measuring_Lights(Current_Protocol)
        program.pre_pulse_light = pre_pulse_light
        program.Ref_Channel = Ref_channel_number(Current_Protocol)
        program.Xe_intensity_value = Xe_intensity_value

        For loopindex As Integer = 0 To Number_Loops(Current_Protocol) - 1 ' VB for loops don't work the same as Java, C#, C, or C++ for loops
            program.Blue_Actinic(loopindex) = M_Blue_actinic(Current_Protocol, loopindex + 1)
            program.Far_Red(loopindex) = M_Far_Red(Current_Protocol, loopindex + 1)
            program.Intensity(loopindex) = Intensity(Current_Protocol, loopindex + 1)
            program.Number_Pulses(loopindex) = Number_Pulses(Current_Protocol, loopindex + 1)
            program.Pre_Delay(loopindex) = Pre_Delay(Current_Protocol, loopindex + 1)
            program.Pre_Pulse(loopindex) = pre_pulse(Current_Protocol, loopindex + 1)
            program.Pre_Pulse_Time(loopindex) = pre_pulse_time(Current_Protocol, loopindex + 1)
            'program.Saturating_Pulse(loopindex) = 
            'program.choke_actinic(loopindex) = m_choke_actinic(Current_Protocol, loopindex + 1)

            program.Take_Data(loopindex) = m_take_data(Current_Protocol, loopindex + 1)
            program.Xe_Flash(loopindex) = Xe_Flash(Current_Protocol, loopindex + 1)
            ' update expected_run_time
            If program.Pre_Delay(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Delay(loopindex))
            If program.Pre_Pulse_Time(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Pulse_Time(loopindex))
            For lightindex As Integer = 0 To Number_Measuring_Lights(Current_Protocol) - 1
                program.Measuring_Interval(lightindex, loopindex) = Measuring_Interval(Current_Protocol, lightindex + 1, loopindex + 1)
                program.Measuring_Light(lightindex) = Measuring_Light(Current_Protocol, lightindex + 1)
                program.choke_actinic(lightindex) = Measuring_Light(Current_Protocol, lightindex + 1)
                program.Primary_Gain(lightindex) = Primary_Gain(Current_Protocol, lightindex + 1)
                program.Reference_Gain(lightindex) = reference_gain(Current_Protocol, lightindex + 1)
                expected_run_time = expected_run_time + _
                                    (Protocol.TimeInSeconds(program.Measuring_Interval(lightindex, loopindex)) _
                                     * program.Number_Pulses(loopindex)) + _
                                    Protocol.TimeInSeconds("20u")
            Next lightindex
        Next loopindex
        Dim data_time_in_doubles(-1) As Double
        Dim programed_successfully As Boolean = cblaster.ProgramProtocol(memory_index, program, data_time_in_doubles)
        If Not programed_successfully Then
            System.Windows.Forms.MessageBox.Show("There was an error while programming protocol " & Current_Protocol & " to the ChrisBlaster (memory address " & memory_index & ").", _
                "Error!")
            System.Windows.Forms.Application.DoEvents()
        End If


    End Sub

    Private Sub Multi_Trace_from_FPGA(ByRef Current_Protocol As Short, ByVal memory_index As Integer, _
                                      ByRef Number_Loops() As Short, ByRef Intensity(,) As Short, ByRef Number_Pulses(,) As Short, _
                                      ByRef Number_Measuring_Lights() As Short, ByRef Measuring_Light(,) As Short, ByRef Choke_Actinic(,) As Short, _
                                      ByRef Measuring_Interval(,,) As String, ByRef Primary_Gain(,) As Short, ByRef reference_gain(,) As Short, _
                                      ByRef Measuring_Pulse_Duration(,) As String, ByRef Gain() As Integer, ByRef In_Channel() As Short, ByRef Ref_Channel_Number() As Short, _
                                        ByRef points As Integer, _
                                      ByRef data_time() As Single, ByRef Xe_Flash(,) As Short, ByRef q_switch(,) As Short, ByVal M_Far_Red(,) As Short, _
                                      ByVal M_Blue_actinic(,) As Short, ByVal pre_pulse(,) As Short, ByVal pre_pulse_time(,) As String, _
                                      ByVal Pre_Delay(,) As String, ByVal m_take_data(,) As Short)

        'Dim ulstat As Short
        'Dim Start_of_Wait, end_of_wait As Single

        ''ulstat = cbDOut(BoardNum, AUXPORT, 0)

        'Trace_Running = True

        ''Dim Time_In_S As Single
        'Dim PP_Line As Short = 0
        ''Dim i, j As Object
        'Dim i As Object
        'Dim ii As Integer
        ''Dim k As Short
        ''Dim AD_Trigger As Short
        'Dim Loop_Return_Index() As Short
        ''Dim Points&
        'Dim CBRate As Integer = 50000
        'Dim CBCount As Integer
        'Dim LowChan, HighChan As Short
        'Dim Options As Integer
        ''Dim ulstat%
        'Dim ss As Short
        'Dim CurCount, CurIndex As Integer
        'Dim Time_Accumulator As Single = 0
        'Dim Measuring_Pulse_Number As Integer = 0
        'Dim Timed_Out As Short
        'Dim Time_Out_Time As Single = 1
        'Dim Each_Measuring_Pulse_Interval() As Single
        'Dim Current_Measuring_Point As Integer = 0
        ''Dim Give_Saturating_Pulse As Short
        'Dim SH_State As Integer
        ''Dim TTT1, TTT2 As Single
        'Dim Raw_Data() As Short

        'ReDim Loop_Return_Index(Number_Loops(Current_Protocol))
        'ReDim Each_Measuring_Pulse_Interval(Number_Measuring_Lights(Current_Protocol))
        ''Dim pbstat As Integer


        'points = 0

        'For i = 1 To Number_Loops(Current_Protocol)
        '    If m_take_data(Current_Protocol, i) > 0 Then
        '        points = points + (Number_Measuring_Lights(Current_Protocol) * Number_Pulses(Current_Protocol, i))
        '    Else
        '        List1.Items.Add("skip data collection found at loop " + Str(i))


        '    End If
        'Next i

        'ReDim data_time(points * 4 + 10)
        'ReDim Data_Volts(points + 10, 4)
        'ReDim Ref_Data_Volts(points * 4 + 10)
        'ReDim Raw_Data(points * 4 + 10) ' dimension an array to hold the input values
        'Start_of_Wait = VB.Timer()
        ''Dim time_temp As String

        '' insert a 1 s delay
        'SH_State = 0 'hold mode

        'ProgressBar1.Value = Primary_Gain(Current_Protocol, 1)
        'sample_gain_label.Text = Primary_Gain(Current_Protocol, 1)
        'ProgressBar2.Value = reference_gain(Current_Protocol, 1)
        'reference_gain_label.Text = reference_gain(Current_Protocol, 1)
        'ProgressBar3.Value = Intensity(Current_Protocol, 1)
        'actinic_label.Text = Intensity(Current_Protocol, 1)


        'Dim expected_run_time As Double = 0 ' This variable keeps track of how long to wait for a protocol to finish (used to detect time-out errors)
        '' Use the ChrisBlaster to run the protocol
        'Dim program As Protocol = New Protocol()

        '' The protocol object holds all of the relevent values of a protocol in one easy package for use by the ChrisBlasterIO object
        '' to handle the nuts & bolts of programming the hardware
        'program.Adc_Gain = Gain
        ''program.Baseline_End = 
        ''program.Baseline_Start = 
        'program.In_Channel = In_Channel
        'program.laser = laser
        'For n_light As Integer = 1 To Number_Measuring_Lights(Current_Protocol)
        '    program.Measuring_Pulse_Duration(n_light) = Measuring_Pulse_Duration(Current_Protocol, n_light)
        'Next n_light

        'program.Number_Loops = Number_Loops(Current_Protocol)
        'program.Number_Measuring_Lights = Number_Measuring_Lights(Current_Protocol)
        'program.pre_pulse_light = pre_pulse_light
        'program.Ref_Channel = use_ref_channel
        'program.Xe_intensity_value = Xe_intensity_value

        'For loopindex As Integer = 0 To Number_Loops(Current_Protocol) - 1 ' VB for loops don't work the same as Java, C#, C, or C++ for loops
        '    program.Blue_Actinic(loopindex) = M_Blue_actinic(Current_Protocol, loopindex + 1)
        '    program.Far_Red(loopindex) = M_Far_Red(Current_Protocol, loopindex + 1)
        '    program.Intensity(loopindex) = Intensity(Current_Protocol, loopindex + 1)
        '    program.Number_Pulses(loopindex) = Number_Pulses(Current_Protocol, loopindex + 1)
        '    program.Pre_Delay(loopindex) = Pre_Delay(Current_Protocol, loopindex + 1)
        '    program.Pre_Pulse(loopindex) = pre_pulse(Current_Protocol, loopindex + 1)
        '    program.Pre_Pulse_Time(loopindex) = pre_pulse_time(Current_Protocol, loopindex + 1)
        '    'program.Saturating_Pulse(loopindex) = 
        '    program.Take_Data(loopindex) = m_take_data(Current_Protocol, loopindex + 1)
        '    program.Xe_Flash(loopindex) = Xe_Flash(Current_Protocol, loopindex + 1)
        '    ' update expected_run_time
        '    If program.Pre_Delay(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Delay(loopindex))
        '    If program.Pre_Pulse_Time(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Pulse_Time(loopindex))
        '    For lightindex As Integer = 0 To Number_Measuring_Lights(Current_Protocol) - 1
        '        program.Measuring_Interval(lightindex, loopindex) = Measuring_Interval(Current_Protocol, lightindex + 1, loopindex + 1)
        '        program.Measuring_Light(lightindex) = Measuring_Light(Current_Protocol, lightindex + 1)
        '        program.choke_actinic(lightindex) = Choke_Actinic(Current_Protocol, lightindex + 1)

        '        program.Primary_Gain(lightindex) = Primary_Gain(Current_Protocol, lightindex + 1)
        '        program.Reference_Gain(lightindex) = reference_gain(Current_Protocol, lightindex + 1)
        '        expected_run_time = expected_run_time + _
        '                            (Protocol.TimeInSeconds(program.Measuring_Interval(lightindex, loopindex)) _
        '                             * program.Number_Pulses(loopindex)) + _
        '                            Protocol.TimeInSeconds("20u")
        '    Next lightindex
        'Next loopindex
        'Dim data_time_in_doubles(-1) As Double
        'data_time_in_doubles = cblaster.GenerateTimeMap(program)
        'ReDim data_time(data_time_in_doubles.Length)
        'For t As Integer = 0 To data_time_in_doubles.Length - 1
        '    data_time(t) = Convert.ToSingle(data_time_in_doubles(t))
        'Next


        'Time_Out_Time = 2 * data_time(data_time.Length - 2) ' time point of last data point


        ''Midnight bug fix:
        'Dim time_now As Double = VB.Timer
        'If (time_now < Start_of_Wait) Then
        '    end_of_wait = Start_of_Wait - time_now - (24 * 60 * 60) ' midnight adjusted code
        'Else
        '    end_of_wait = Start_of_Wait - VB.Timer() ' original code
        'End If

        ''List1.AddItem "points:" & Str$(Points&)

        ''************************

        '' Collect the values with cbAInScan%()
        '' Parameters:
        ''   BoardNum%   :the number used by CB.CFG to describe this board
        ''   LowChan%    :the first channel of the scan
        ''   HighChan%   :the last channel of the scan
        ''   CBCount&    :the total number of A/D samples to collect
        ''   CBRate&     :sample rate
        ''   Gain        :the gain for the board
        ''   ADData%     :the array for the collected data values
        ''   Options     :data collection options


        'Time_Out_Time = Time_Out_Time * 5


        ''CBRate = 3 ' sampling rate (samples per second)

        'LowChan = 0 'In_Channel '0                   ' first channel to acquire
        'HighChan = 3 'In_Channel ' 0                  ' last channel to acquire
        'number_of_channels_to_scan = HighChan - LowChan ' used later in data analyses

        'Const FirstPoint As Integer = 0 ' set first element in buffer to transfer to array

        '' total number of data points to collect

        'CBCount = points * number_of_channels_to_scan  'highchannel-lowchannel
        ''If points > 10000 Then

        'MemHandle = cbWinBufAlloc(points * number_of_channels_to_scan + 10) ' set aside memory to hold data
        ''End If
        ''      MemHandle = cbWinBufAlloc(CBCount) ' set aside memory to hold data
        'ulstat = cbSetTrigger(BoardNum, TRIGPOSEDGE, 0, 0.5)  'changes high from 1 to 2

        ''Options = EXTTRIGGER + BACKGROUND + NOCONVERTDATA '+ SINGLEIO
        'Options = EXTCLOCK + BACKGROUND + NOCONVERTDATA '+ SINGLEIO

        '' return data as 16-bit values
        '' collect data in BACKGROUND mode
        '' use NOCONVERTDATA if using 16 bit board


        ''        List1.Items.Add(" trace started")
        'System.Windows.Forms.Application.DoEvents()

        'If MemHandle = 0 Then Stop ' check that a handle to a memory buffer exists

        'ulstat = cbAInScan(BoardNum, LowChan, HighChan, CBCount, CBRate, Gain, MemHandle, Options)


        ''Call Hold_on_There(3)

        ''Call outt(BaseAddress + 6, 0)  'not needed for this version of PB?
        ''Call outt(BaseAddress + 7, 0) '7)


        ''put the device in run mode
        ''Call outt(BaseAddress + 1, 0) '7)
        'If (Not cblaster.StartRunningProtocol(memory_index)) Then
        '    System.Windows.Forms.MessageBox.Show("Protocol " & Current_Protocol & " from memory index " & memory_index & " failed to run or timed out.")
        '    System.Windows.Forms.Application.DoEvents()
        'End If


        'ulstat = cbGetStatus(BoardNum, ss, CurCount, CurIndex, AIFUNCTION)

        'Timed_Out = False
        'Time_Out_Time = VB.Timer() + Time_Out_Time
        ''Midnight bug fix:
        '' First, find the correct ending wait time
        'Dim midnight_mark As Boolean = False
        'If (Time_Out_Time > (24 * 60 * 60)) Then
        '    Time_Out_Time = Time_Out_Time - (24 * 60 * 60)
        '    midnight_mark = True ' tell following code not to trigger "timed out" until after midnight
        'End If

        'Label13.Text = CurCount & " / " & CBCount

        'While CurCount < CBCount And Timed_Out = False
        '    If (VB.Timer() > Time_Out_Time And midnight_mark = False) Then Timed_Out = True ' don't time-out before midnight
        '    If (midnight_mark = True And VB.Timer() < 60) Then midnight_mark = False ' indicate that midnight has passed

        '    System.Windows.Forms.Application.DoEvents()
        '    ulstat = cbGetStatus(BoardNum, ss, CurCount, CurIndex, AIFUNCTION)
        '    Label13.Text = CurCount & " / " & CBCount
        '    'System.Windows.Forms.Application.DoEvents()
        'End While

        'Label3.Text = CurCount & " * " & CBCount
        'System.Windows.Forms.Application.DoEvents()
        ''List1.Items.Add(CBCount & "CC:" & CurCount)
        ''ulstat = cbDOut(BoardNum, AUXPORT, &HFFS)

        'If Timed_Out = True Then
        '    ulstat = cbStopBackground(BoardNum, AIFUNCTION)
        '    System.Windows.Forms.MessageBox.Show("Timed out!")
        '    System.Windows.Forms.Application.DoEvents()

        'Else
        '    List1.Items.Add(" trace done")
        '    'Scrolling text box fix: make the box jump to the bottom after adding an item
        '    List1.TopIndex = List1.Items.Count - 1
        'End If


        'ulstat = cbStopBackground(BoardNum, AIFUNCTION) 'clears variables and flags after normal termination

        ''Transfer the data from the memory buffer set up by Windows to an array for use by Visual Basic

        'ulstat = cbWinBufToArray(MemHandle, Raw_Data(0), FirstPoint, CBCount)

        ''If ulstat% <> 0 Then Stop
        '' convert RAW data into VOLTS

        ''        For i = 0 To (points - 1) * 2 Step 2
        '' ulstat = cbToEngUnits(BoardNum, Gain, Raw_Data(i), TTT1)
        ''ulstat = cbToEngUnits(BoardNum, Gain, Raw_Data(i + 1), TTT2)
        ''Data_Volts(i / 2) = (TTT1 / TTT2)
        ''Next i
        '' Data_Volts is zero-indexed
        '' Raw_Data is zero_indexed

        'For i = 0 To (points - 1)
        '    For ii = 0 To number_of_channels_to_scan - 1 'this change may allow us to scan more than 4 channels
        '        ulstat = cbToEngUnits(BoardNum, Gain, Raw_Data((i * 4) + ii), Data_Volts(i, ii))
        '    Next ii
        'Next i


        'Trace_Running = False
    End Sub

    Private Sub Multi_Trace(ByRef Current_Protocol As Short, ByRef Number_Loops() As Short, ByRef Intensity(,) _
                            As Short, ByRef Number_Pulses(,) As Short, ByRef Number_Measuring_Lights() As Short, _
                            ByRef Measuring_Light(,) As Short, ByRef choke_actinic(,) As Short, ByRef Measuring_Interval(,,) As String, _
                            ByRef Primary_Gain(,) As Short, ByRef reference_gain(,) As Short, ByRef Measuring_Pulse_Duration(,) As String, _
                            ByRef Gain() As Integer, ByRef In_Channel() As Short, ByRef Ref_channel_number() As Short, ByRef points As Integer, ByRef data_time() As Single, _
                            ByRef Xe_Flash(,) As Short, ByRef q_switch(,) As Short, ByVal M_Far_Red(,) As Short, ByVal M_Blue_actinic(,) _
                            As Short, ByVal pre_pulse(,) As Short, ByVal pre_pulse_time(,) As String, ByVal Pre_Delay(,) As String, _
                            ByVal m_take_data(,) As Short)

        halt.Enabled = True

        ' Dim reenable_timers As Boolean = False

        ' reenable_timers = False

        ' If Form3.Timer1.Enabled = True And Form3.Timer1.Enabled = True Then
        '     reenable_timers = True
        '     Form3.Timer1.Enabled = False ' disable the monochromator timers

        '     Form3.Timer2.Enabled = False ' disable the monochromator timers
        ' End If



        Dim ulstat As Short
        Dim Start_of_Wait, end_of_wait As Single

        'ulstat = cbDOut(BoardNum, AUXPORT, 0)

        Trace_Running = True

        'Dim Time_In_S As Single
        Dim PP_Line As Short = 0
        'Dim i, j As Object
        Dim i As Object
        Dim ii As Integer
        'Dim k As Short
        'Dim AD_Trigger As Short
        Dim Loop_Return_Index() As Short
        'Dim Points&
        Dim CBRate As Integer = 50000
        Dim CBCount As Integer
        Dim LowChan, HighChan As Short
        Dim Options As Integer
        'Dim ulstat%
        Dim ss As Short
        Dim CurCount, CurIndex As Integer
        Dim Time_Accumulator As Single = 0
        Dim Measuring_Pulse_Number As Integer = 0
        Dim Timed_Out As Short
        Dim Time_Out_Time As Single = 1
        Dim Each_Measuring_Pulse_Interval() As Single
        Dim Current_Measuring_Point As Integer = 0
        'Dim Give_Saturating_Pulse As Short
        Dim SH_State As Integer
        'Dim TTT1, TTT2 As Single
        Dim Raw_Data() As Short

        ReDim Loop_Return_Index(Number_Loops(Current_Protocol))
        ReDim Each_Measuring_Pulse_Interval(Number_Measuring_Lights(Current_Protocol))
        'Dim pbstat As Integer

        Dim tt As Integer = m_take_data.Length()

        points = 0

        For i = 1 To Number_Loops(Current_Protocol)
            If m_take_data(Current_Protocol, i) > 0 Then
                points = points + (Number_Measuring_Lights(Current_Protocol) * Number_Pulses(Current_Protocol, i))
            End If
        Next i

        'List1.Items.Add("protocol: " + Str(Current_Protocol) + " points = " + Str(points))
        LowChan = 0 'In_Channel '0                   ' first channel to acquire
        HighChan = 3 'In_Channel ' 0                  ' last channel to acquire
        number_of_channels_to_scan = HighChan - LowChan + 1 ' used later in data analyses

        ReDim data_time(points * number_of_channels_to_scan + 10)
        ReDim Data_Volts(points + 10, number_of_channels_to_scan)
        ReDim Ref_Data_Volts(points * number_of_channels_to_scan + 10)
        ReDim Raw_Data(points * number_of_channels_to_scan + 10) ' dimension an array to hold the input values

        Start_of_Wait = VB.Timer()
        'Dim time_temp As String

        ' insert a 1 s delay
        SH_State = 0 'hold mode

        ProgressBar1.Value = Primary_Gain(Current_Protocol, 1)
        sample_gain_label.Text = Primary_Gain(Current_Protocol, 1)
        ProgressBar2.Value = reference_gain(Current_Protocol, 1)
        reference_gain_label.Text = reference_gain(Current_Protocol, 1)
        ProgressBar3.Value = Intensity(Current_Protocol, 1)
        actinic_label.Text = Intensity(Current_Protocol, 1)
        Dim ix As Short
        'Jabber_jabber.gain_data_list.Items.Add("after: ")
        'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
        '    Jabber_jabber.gain_data_list.Items.Add(m_Primary_Gain(Current_Protocol, ix))
        'Next
        'For ix = 1 To m_Number_Measuring_Lights(Current_Protocol)
        '    Jabber_jabber.gain_data_list.Items.Add(m_reference_gain(Current_Protocol, ix))
        'Next

        Dim expected_run_time As Double = 0 ' This variable keeps track of how long to wait for a protocol to finish (used to detect time-out errors)
        ' Use the ChrisBlaster to run the protocol
        Dim program As Protocol = New Protocol()

        ' The protocol object holds all of the relevent values of a protocol in one easy package for use by the ChrisBlasterIO object
        ' to handle the nuts & bolts of programming the hardware
        program.Adc_Gain = Gain(Current_Protocol)
        'program.Baseline_End = 
        'program.Baseline_Start = 
        program.In_Channel = In_Channel(Current_Protocol)
        program.Ref_Channel = Ref_channel_number(Current_Protocol)

        program.laser = laser
        program.Number_Measuring_Lights = Number_Measuring_Lights(Current_Protocol)

        For n_light As Integer = 1 To Number_Measuring_Lights(Current_Protocol)
            'PROBLEM idemntified:
            ' when using nLIGHTS>5 crashes. Most likely caused by FUd 1 versus 0 indexing
            ' Note: I figured this out once, and it seems to be in the Measuring_Pulse_Duration 
            'DMK_update: program.Measuring_Pulse_Duration(n_light - 1) = Measuring_Pulse_Duration(Current_Protocol, n_light) 
            'a possible problem! changed program.Measuring_Pulse_Duration(n_light) to program.Measuring_Pulse_Duration(n_light - 1)
            'Note: program.Measuring_Pulse_Duration is zero-indexed, Measuring_Pulse_Duration(Current_Protocol, n_light) is 1-indexed!!
            '                   I thus changed to offset the index when transferring the values into the 
            '                   object by subtracting one from the index program.Measuring_Pulse_Duration(n_light - 1)
            program.Measuring_Pulse_Duration(n_light - 1) = Measuring_Pulse_Duration(Current_Protocol, n_light) 'a possible problem! changed program.Measuring_Pulse_Duration(n_light) to program.Measuring_Pulse_Duration(n_light - 1)
        Next
        program.Number_Loops = Number_Loops(Current_Protocol)

        program.pre_pulse_light = pre_pulse_light
        program.Xe_intensity_value = Xe_intensity_value
        '@@
        For loopindex As Integer = 0 To Number_Loops(Current_Protocol) - 1 ' VB for loops don't work the same as Java, C#, C, or C++ for loops
            program.Blue_Actinic(loopindex) = M_Blue_actinic(Current_Protocol, loopindex + 1)

            program.Far_Red(loopindex) = M_Far_Red(Current_Protocol, loopindex + 1)
            program.Intensity(loopindex) = Intensity(Current_Protocol, loopindex + 1)
            program.Number_Pulses(loopindex) = Number_Pulses(Current_Protocol, loopindex + 1)
            program.Pre_Delay(loopindex) = Pre_Delay(Current_Protocol, loopindex + 1)
            program.Pre_Pulse(loopindex) = pre_pulse(Current_Protocol, loopindex + 1)
            program.Pre_Pulse_Time(loopindex) = pre_pulse_time(Current_Protocol, loopindex + 1)
            'program.Saturating_Pulse(loopindex) = 
            program.Take_Data(loopindex) = m_take_data(Current_Protocol, loopindex + 1)
            program.Xe_Flash(loopindex) = Xe_Flash(Current_Protocol, loopindex + 1)
            ' update expected_run_time
            If program.Pre_Delay(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Delay(loopindex))
            If program.Pre_Pulse_Time(loopindex) <> "" Then expected_run_time = expected_run_time + Protocol.TimeInSeconds(program.Pre_Pulse_Time(loopindex))
            For lightindex As Integer = 0 To Number_Measuring_Lights(Current_Protocol) - 1

                program.Measuring_Interval(lightindex, loopindex) = Measuring_Interval(Current_Protocol, lightindex + 1, loopindex + 1)

                'program.m_wasp(lightindex) = M_Wasp(Current_Protocol, lightindex + 1)

                program.Measuring_Light(lightindex) = Measuring_Light(Current_Protocol, lightindex + 1)

                program.choke_actinic(lightindex) = choke_actinic(Current_Protocol, lightindex + 1)


                program.Primary_Gain(lightindex) = Primary_Gain(Current_Protocol, lightindex + 1)
                'gain_data_list.Items.Add(" !!pg: " + Str(Primary_Gain(Current_Protocol, lightindex + 1)))

                program.Reference_Gain(lightindex) = reference_gain(Current_Protocol, lightindex + 1)

                expected_run_time = expected_run_time + _
                                    (Protocol.TimeInSeconds(program.Measuring_Interval(lightindex, loopindex)) _
                                     * program.Number_Pulses(loopindex)) + _
                                    Protocol.TimeInSeconds("20u")
                'Jabber_jabber.List1.Items.Add("protocol: " + Str(Current_Protocol) + "%%%: " + program.Measuring_Pulse_Duration(lightindex)) ' + " " + Str(program.Measuring_Light(lightindex)) + "  " + Str(program.choke_actinic(lightindex)) + " " + Str(program.Primary_Gain(lightindex)) + " " + Str(program.Reference_Gain(lightindex)))

            Next lightindex
        Next loopindex
        'Jabber_jabber.List1.Items.Add("expected run time= " + Str(expected_run_time))
        Dim data_time_in_doubles(-1) As Double
        '*****************************************this is where the cblaster is programmed************************************

        Dim programed_successfully As Boolean = cblaster.ProgramProtocol(Current_Protocol, program, data_time_in_doubles)

        '*****************************************above is where the cblaster is programmed************************************

        If programed_successfully Then
            ReDim data_time(data_time_in_doubles.Length)
            For t As Integer = 0 To data_time_in_doubles.Length - 1
                data_time(t) = Convert.ToSingle(data_time_in_doubles(t))
            Next
        End If
        Jabber_jabber.List1.Items.Add(" trace started")

        Time_Out_Time = data_time(data_time.Length - 2) + 1 ' time point of last data point and multiple by 2 for good measure

        Jabber_jabber.List1.Items.Add("Time_Out_Time = " + Str(Time_Out_Time))


        'Midnight bug fix:
        Dim time_now As Double = VB.Timer
        If (time_now < Start_of_Wait) Then
            end_of_wait = Start_of_Wait - time_now - (24 * 60 * 60) ' midnight adjusted code
        Else
            end_of_wait = Start_of_Wait - VB.Timer() ' original code
        End If


        'Call outt(BaseAddress% + 6, 0)
        'Call outt(BaseAddress% + 7, 7)


        'List1.AddItem "points:" & Str$(Points&)

        '************************

        ' Collect the values with cbAInScan%()
        ' Parameters:
        '   BoardNum%   :the number used by CB.CFG to describe this board
        '   LowChan%    :the first channel of the scan
        '   HighChan%   :the last channel of the scan
        '   CBCount&    :the total number of A/D samples to collect
        '   CBRate&     :sample rate
        '   Gain        :the gain for the board
        '   ADData%     :the array for the collected data values
        '   Options     :data collection options


        Time_Out_Time = Time_Out_Time * 2


        'CBRate = 3 ' sampling rate (samples per second)

        'LowChan = 0 'In_Channel '0                   ' first channel to acquire
        'HighChan = 3 'In_Channel ' 0                  ' last channel to acquire

        Const FirstPoint As Integer = 0 ' set first element in buffer to transfer to array

        ' total number of data points to collect

        CBCount = points * number_of_channels_to_scan  'highchannel-lowchannel
        'If points > 10000 Then

        MemHandle = cbWinBufAlloc(points * number_of_channels_to_scan + 10) ' set aside memory to hold data
        'End If
        '      MemHandle = cbWinBufAlloc(CBCount) ' set aside memory to hold data
        ulstat = cbSetTrigger(BoardNum, TRIGPOSEDGE, 0, 0.5)  'changes high from 1 to 2

        'Options = EXTTRIGGER + BACKGROUND + NOCONVERTDATA '+ SINGLEIO
        Options = EXTCLOCK + BACKGROUND + NOCONVERTDATA '+ SINGLEIO

        ' return data as 16-bit values
        ' collect data in BACKGROUND mode
        ' use NOCONVERTDATA if using 16 bit board


        '        List1.Items.Add(" trace started")
        System.Windows.Forms.Application.DoEvents()

        If MemHandle = 0 Then
            Jabber_jabber.List1.Items.Add("Memhandle failed")
            Stop ' check that a handle to a memory buffer exists
        End If

        ulstat = cbAInScan(BoardNum, LowChan, HighChan, CBCount, CBRate, Gain(Current_Protocol), MemHandle, Options)


        'Call Hold_on_There(3)

        'Call outt(BaseAddress + 6, 0)  'not needed for this version of PB?
        'Call outt(BaseAddress + 7, 0) '7)


        'put the device in run mode
        'Call outt(BaseAddress + 1, 0) '7)
        If (Not cblaster.StartRunningProtocol(Current_Protocol)) Then
            System.Windows.Forms.MessageBox.Show("Protocol " & Current_Protocol & " failed to run or timed out.")
            System.Windows.Forms.Application.DoEvents()
        End If


        ulstat = cbGetStatus(BoardNum, ss, CurCount, CurIndex, AIFUNCTION)

        Timed_Out = False
        Time_Out_Time = VB.Timer() + Time_Out_Time
        'Midnight bug fix:
        ' First, find the correct ending wait time
        Dim midnight_mark As Boolean = False
        If (Time_Out_Time > (24 * 60 * 60)) Then
            Time_Out_Time = Time_Out_Time - (24 * 60 * 60)
            midnight_mark = True ' tell following code not to trigger "timed out" until after midnight
        End If

        'List1.Items.Add(CurCount & " / " & CBCount)
        'Dim count As Integer
        'Dim ulstat2
        If spec_mode = True Then
            'ulstat = cbCLoad(1, 1, 0)
        End If

        While CurCount < CBCount And Timed_Out = False
            If (VB.Timer() > Time_Out_Time And midnight_mark = False) Then
                Timed_Out = True ' don't time-out before midnight
                Jabber_jabber.List1.Items.Add("timed out (midnight=False)")
            End If
            If (midnight_mark = True And VB.Timer() < 60) Then
                Jabber_jabber.List1.Items.Add("timed out at not mid")
                midnight_mark = False ' indicate that midnight has passed

            End If
            System.Windows.Forms.Application.DoEvents()
            ulstat = cbGetStatus(BoardNum, ss, CurCount, CurIndex, AIFUNCTION)
            Label13.Text = CurCount & " / " & CBCount


            'Label4.Text = count
            If spec_mode = True Then
                'ulstat2 = cbCIn(1, 1, count)
                'Set_Lambda(count * delta_lambda + start_lambda)
                'Label4.Text = delta_lambda
                'Set_Lambda(start_lambda)
            End If
            System.Windows.Forms.Application.DoEvents()
        End While

        Label13.Text = CurCount & " ! " & CBCount
        System.Windows.Forms.Application.DoEvents()
        'List1.Items.Add(CBCount & "CC:" & CurCount)
        'ulstat = cbDOut(BoardNum, AUXPORT, &HFFS)

        If Timed_Out = True Then
            ulstat = cbStopBackground(BoardNum, AIFUNCTION)
            System.Windows.Forms.MessageBox.Show("Timed out!")
            System.Windows.Forms.Application.DoEvents()

        Else
            Jabber_jabber.List1.Items.Add(" trace done")
            'List1.Items.Add(Str(CurCount) + " * " + Str(CBCount))
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            Jabber_jabber.List1.TopIndex = List1.Items.Count - 1
        End If


        ulstat = cbStopBackground(BoardNum, AIFUNCTION) 'clears variables and flags after normal termination

        'Transfer the data from the memory buffer set up by Windows to an array for use by Visual Basic

        ulstat = cbWinBufToArray(MemHandle, Raw_Data(0), FirstPoint, CBCount)

        'If ulstat% <> 0 Then Stop
        ' convert RAW data into VOLTS

        '        For i = 0 To (points - 1) * 2 Step 2
        ' ulstat = cbToEngUnits(BoardNum, Gain, Raw_Data(i), TTT1)
        'ulstat = cbToEngUnits(BoardNum, Gain, Raw_Data(i + 1), TTT2)
        'Data_Volts(i / 2) = (TTT1 / TTT2)
        'Next i

        For i = 0 To (points - 1)
            For ii = 0 To number_of_channels_to_scan - 1
                ulstat = cbToEngUnits(BoardNum, Gain(Current_Protocol), Raw_Data((i * 4) + ii), Data_Volts(i, ii))
                'List1.Items.Add("i:" & Str(i) & " ii:" & Str(ii) & " rdi:" & Str((i * 4) + ii) & "d:" & Str(Data_Volts(i, ii)))
            Next ii
        Next i

        Dim testvar As Integer = Raw_Data.Length()
        Dim testvar2 As Integer = Data_Volts.Length()

        Trace_Running = False


        ' If reenable_timers = True Then
        '     Form3.Timer1.Enabled = True ' disable the monochromator timers

        '     Form3.Timer2.Enabled = True ' disable the monochromator timers
        ' End If


    End Sub



    Public Sub dmk_exit_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles dmk_exit.Click
        End
    End Sub

    Private Sub DMK_v1_Load(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles MyBase.Load
        Dim ulstat As Short
        'Dim i, ii As Short
        Dim i As Short
        'Dim PBstat As Short
        'current_lambda = My.Settings.lambda
        LambdaText.Text = current_lambda

        BaseAddress = 0

        '  declare revision level of Universal Library

        ulstat = cbDeclareRevision(CURRENTREVNUM)

        ' Initiate error handling
        '  activating error handling will trap errors like
        '  bad channel numbers and non-configured conditions.
        '  Parameters:
        '    PRINTALL    :all warnings and errors encountered will be printed
        '    DONTSTOP    :if an error is encountered, the program will not stop,
        '                  errors must be handled locally

        ulstat = cbErrHandling(PRINTALL, STOPALL)
        If ulstat <> 0 Then Stop

        ' If cbErrHandling% is set for STOPALL or STOPFATAL during the program
        ' design stage, Visual Basic will be unloaded when an error is encountered.
        ' We suggest trapping errors locally until the program is ready for compiling
        ' to avoid losing unsaved data during program design.  This can be done by
        ' setting cbErrHandling options as above and checking the value of ULStat%
        ' after a call to the library. If it is not equal to 0, an error has occurred.

        'FileOpen(1, "c:\testx.dat", OpenMode.Output)

        'For i = 1 To 100
        ' ii = 1
        ' Print(1, i & Chr(9) & System.Math.Sin(i * ii / 100) & Chr(9) & System.Math.Cos(i * ii / 100))
        ' For ii = 2 To 5
        ' Print(1, Chr(9) & i & Chr(9) & System.Math.Sin(i * ii / 100) & Chr(9) & System.Math.Cos(i * ii / 100))
        ' Next ii
        ' PrintLine(1, "")
        ' Next i
        ' FileClose(1)
        Dim j_m_plot_trace(3, 3) As Short
        j_m_plot_trace(0, 1) = 1

        For i = 1 To 6
            'Jabber_jabber.List1.Items.Add(Str(j_m_plot_trace(0, 1)))
            Call Plot_File(rootdata_directory & "test.dat", (i), (i And 1), False, 0, 0, j_m_plot_trace)
        Next i


        FileOpen(1, rootdata_directory & "temp.dat", OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        PrintLine(1, "")
        FileClose(1)
        'Threading.Thread.sleep(10)

        '        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)
        '        D_Out_Mask = 0
        '        ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)

        stir_control(0)
    End Sub


    Sub Hold_on_There(ByRef Wait_Time As Single)
        Dim Start_of_Wait As Single

        Start_of_Wait = VB.Timer()
        Bug_Out = False

        'Midnight bug fix:
        ' Adjust the waiting time by subtracting a day if the wait should end tomorrow
        If ((Wait_Time + Start_of_Wait) > (24 * 60 * 60)) Then
            Wait_Time = Wait_Time - (24 * 60 * 60)
        End If

        While VB.Timer() - Start_of_Wait < Wait_Time And Bug_Out = False And Halt_Script = False

            Label2.Text = CStr(Int(10 * (Wait_Time - (VB.Timer() - Start_of_Wait))) / 10)

            System.Windows.Forms.Application.DoEvents()

        End While

    End Sub
    Sub combineJSON(ByVal Base_File_Name As String)
        Dim j As String
        Jabber_jabber.List1.Items.Clear()
        Dim bfolder() As String = Base_File_Name.Split("\")
        Dim endtxt = bfolder(bfolder.Length - 1)

        Dim foldd() As String = Base_File_Name.Split(endtxt(0))
        Dim useThisFolder = foldd(0)
        Dim useThisBaseFile = foldd(1)
        FileOpen(1, Base_File_Name + "combined.jsonx", OpenMode.Append, OpenAccess.Default, OpenShare.Shared)
        Print(1, "[")
        Dim combineFiles = My.Computer.FileSystem.GetFiles(useThisFolder, FileIO.SearchOption.SearchTopLevelOnly, endtxt + "*.json")
        Dim numFiles As Short = combineFiles.Count

        'For Each foundFile As String In combineFiles
        For index As Short = 0 To numFiles - 1
            Dim foundFile As String = combineFiles.item(index)
            'Print(1, "{""file_name"":""")
            bfolder = (foundFile.Split("\"))
            Dim bendtxt As String = bfolder(bfolder.Length - 1)
            'Print(1, bendtxt + """,")
            Jabber_jabber.gain_data_list.Items.Add(foundFile)
            'FileClose(22)
            FileOpen(22, foundFile, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
            While (EOF(22) = False)
                j = LineInput(22)
                Print(1, j)
            End While
            FileClose(22)
            My.Computer.FileSystem.DeleteFile(foundFile)
            If index < numFiles - 1 Then
                Print(1, ",")
            End If
        Next
        Print(1, "]")
        FileClose(1)
        My.Computer.FileSystem.RenameFile(Base_File_Name + "combined.jsonx", endtxt + "combined.json")

    End Sub
    Sub saveJSON(ByRef I_File_Name As String, ByRef Aux_File_Name As String, ByRef Current_Protocol As Short, ByRef Number_Loops() As Short, ByRef Intensity(,) _
                            As Short, ByRef Number_Pulses(,) As Short, ByRef Number_Measuring_Lights() As Short, _
                            ByRef Measuring_Light(,) As Short, ByRef choke_actinic(,) As Short, ByRef Measuring_Interval(,,) As String, _
                            ByRef Primary_Gain(,) As Short, ByRef reference_gain(,) As Short, ByRef Measuring_Pulse_Duration(,) As String, _
                            ByRef Gain() As Integer, ByRef In_Channel() As Short, ByRef Ref_channel_number() As Short, _
                            ByRef Xe_Flash(,) As Short, ByRef q_switch(,) As Short, ByVal M_Far_Red(,) As Short, ByVal M_Blue_actinic(,) As Short, _
                            ByRef trace_note() As String, ByRef trace_label() As String, ByRef M_Measuring_Light_Names(,) As String, _
                            ByRef time_mode() As String, ByRef Start_Time As Single, ByRef current_trace As Short, ByVal protocol_label() As String, ByRef current_lambda As Single)

        Dim Number_columns As Object
        'Dim Number_Measuring_Lights As Short
        Dim i As Integer
        Dim j As String
        Dim p, pp As Integer
        'Dim ppp, p, pp, light As Object
        'Dim Column_N As Short
        Dim index2 As Integer
        Dim factor As Single

        Dim Delta_Data_X(,) As Object
        Dim Delta_Data_Y(,) As Object
        Dim Delta_Data_ref(,) As Object
        Dim Delta_Data_A(,) As Object

        Dim Aux_Delta_Data_X(,) As Object
        Dim Aux_Delta_Data_Y(,) As Object
        Dim Aux_Delta_Data_ref(,) As Object
        Dim Aux_Delta_Data_A(,) As Object


        Dim temp As Single = 0
        'Dim I0_Y As Single
        'Dim I0_Data_Y() As Single
        'Dim Aux_I0_Data_Y() As Single

        ' figure out how many columns are in I_FILE
        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        Dim tttt() = j.Split(Chr(9))

        FileClose(1)

        Number_columns = Number_Measuring_Lights(Current_Protocol)


        'FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        'i = 0

        'While EOF(1) = False
        '    Input(1, temp)
        '    Input(1, temp)
        '    Input(1, temp)
        '    Input(1, temp)
        '    i = i + 1
        'End While
        'FileClose(1)
        'Threading.Thread.sleep(10)

        'dimension the temporary arrays
        Dim number_pts As Integer = 0
        Dim iii As Integer
        For iii = 1 To Number_Loops(Current_Protocol)
            number_pts = number_pts + Number_Pulses(Current_Protocol, iii) '* Number_Measuring_Lights(Current_Protocol)


        Next

        ReDim Delta_Data_X(Number_columns, number_pts)
        ReDim Delta_Data_Y(Number_columns, number_pts)
        'ReDim I0_Data_Y(Number_columns, number_pts)
        ReDim Delta_Data_ref(Number_columns, number_pts)
        ReDim Delta_Data_A(Number_columns, number_pts)


        ' input the RAW data from I_FILE

        Dim ix As Integer = 0

        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        For i = 0 To number_pts - 1
            j = LineInput(1)
            tttt = j.Split(Chr(9))

            For iiii As Integer = 0 To Number_columns - 1
                Dim ii As Short = iiii * 4
                Delta_Data_X(iiii, i) = CSng(tttt(ii))
                Delta_Data_Y(iiii, i) = CSng(tttt(ii + 1))
                Delta_Data_ref(iiii, i) = CSng(tttt(ii + 2))
                Delta_Data_A(iiii, i) = CSng(tttt(ii + 3))

                'Input(1, Delta_Data_X(i))
                'Input(1, Delta_Data_Y(i))
                'Input(1, Delta_Data_ref(i))
                'Input(1, Delta_Data_A(i))

            Next
        Next


        FileClose(1)


        'FileOpen(1, Aux_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        'i = 0

        'While EOF(1) = False
        '    Input(1, temp)
        '    Input(1, temp)
        '    Input(1, temp)
        '    Input(1, temp)
        '    i = i + 1
        'End While
        'FileClose(1)
        'Threading.Thread.sleep(10)

        'dimension the temporary arrays

        ReDim Aux_Delta_Data_X(Number_columns, number_pts)
        ReDim Aux_Delta_Data_Y(Number_columns, number_pts)
        ReDim Aux_Delta_Data_ref(Number_columns, number_pts)
        ReDim Aux_Delta_Data_A(Number_columns, number_pts)


        ' input the RAW data from I_FILE

        i = 0

        FileOpen(1, Aux_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        'LineInput()
        'Dim iii As Short
        For i = 0 To number_pts - 1
            j = LineInput(1)
            tttt = j.Split(Chr(9))

            For iiii As Integer = 0 To Number_columns - 1
                Dim ii As Short = iiii * 4
                Aux_Delta_Data_X(iiii, i) = CSng(tttt(ii))
                Aux_Delta_Data_Y(iiii, i) = CSng(tttt(ii + 1))
                Aux_Delta_Data_ref(iiii, i) = CSng(tttt(ii + 2))
                Aux_Delta_Data_A(iiii, i) = CSng(tttt(ii + 3))

                'Input(1, Delta_Data_X(i))
                'Input(1, Delta_Data_Y(i))
                'Input(1, Delta_Data_ref(i))
                'Input(1, Delta_Data_A(i))

            Next
        Next
        'For i = 1 To number_pts
        '    For ii As Integer = 1 To Number_columns
        '        Input(1, Aux_Delta_Data_X(i))
        '        Input(1, Aux_Delta_Data_Y(i))
        '        Input(1, Aux_Delta_Data_ref(i))
        '        Input(1, Aux_Delta_Data_A(i))
        '        i = i + 1

        '    Next
        'Next
        'Input(1, Aux_Delta_Data_X(i))
        'Input(1, Aux_Delta_Data_Y(i))
        'Input(1, Aux_Delta_Data_ref(i))
        'Input(1, Aux_Delta_Data_A(i))
        'i = i + 1
        'End While

        FileClose(1)


        FileOpen(1, I_File_Name + ".json", OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        Print(1, "{""trace_label"":""")
        Print(1, trace_label(current_trace))
        Print(1, """,")

        Print(1, """trace_number"":")
        Print(1, current_trace)
        Print(1, ",")

        Print(1, """start_time"":")
        Print(1, Start_Time)
        Print(1, ",")

        Print(1, """note"":""")
        Print(1, trace_note(Current_Protocol))
        Print(1, """,")

        Print(1, """wl"":")
        Print(1, current_lambda)
        Print(1, ",")


        'Print(1, """data"":{")

        For pp = 0 To Number_columns - 1
            'Print(1, """" & M_Measuring_Light_Names(Current_Protocol, pp) & """:{")
            Dim n As String = M_Measuring_Light_Names(Current_Protocol, pp)
            'Print(1, """time"":[")
            Print(1, """" + n + "_time"":[")
            For p = 0 To number_pts - 2
                Print(1, Delta_Data_X(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Delta_Data_X(pp, p) & "],")
            Print(1, """" + n + "_I"":[")
            For p = 0 To number_pts - 2
                Print(1, Delta_Data_Y(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Delta_Data_Y(pp, p) & "],")

            'Print(1, """I0"":[")
            Print(1, """" + n + "_I0"":[")
            For p = 0 To number_pts - 2
                Print(1, Delta_Data_ref(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Delta_Data_ref(pp, p) & "],")

            'Print(1, """calc"":[")
            Print(1, """" + n + "_calc"":[")
            For p = 0 To number_pts - 2
                Print(1, Delta_Data_A(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Delta_Data_A(pp, p) & "],")

            'Next pp

            '' now do the aux data
            ''Print(1, """aux_data"":{")
            'For pp = 0 To Number_columns - 1
            '    Print(1, """" & M_Measuring_Light_Names(Current_Protocol, pp) & """:{")

            'Print(1, """auxtime"":[")
            Print(1, """" + n + "auxtime"":[")
            For p = 0 To number_pts - 2
                Print(1, Aux_Delta_Data_X(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Aux_Delta_Data_X(pp, p) & "],")
            Print(1, """" + n + "_Ia"":[")
            'Print(1, """auxI"":[")
            For p = 0 To number_pts - 2
                Print(1, Aux_Delta_Data_Y(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Aux_Delta_Data_Y(pp, p) & "],")

            'Print(1, """auxI0"":[")
            Print(1, """" + n + "_I0a"":[")
            For p = 0 To number_pts - 2
                Print(1, Aux_Delta_Data_ref(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Aux_Delta_Data_ref(pp, p) & "],")

            'Print(1, """auxcalc"":[")
            Print(1, """" + n + "_calca"":[")
            For p = 0 To number_pts - 2
                Print(1, Aux_Delta_Data_A(pp, p) & ",") ' & Chr(9) & Delta_Data_Y(p) & Chr(9) & Delta_Data_ref(p) & Chr(9) & Delta_Data_A(p))
            Next p
            Print(1, Aux_Delta_Data_A(pp, p) & "],")


        Next pp


        'Print(1, """protocol"":{")

        Print(1, """protocol_number"":" & Current_Protocol & ",")

        Print(1, """protocol_label"":""" & protocol_label(Current_Protocol) & """,")

        Print(1, """number_loops"":" & Number_Loops(Current_Protocol) & ",")

        Print(1, """sample_channel"":" & In_Channel(Current_Protocol) & ",")

        Print(1, """reference_channel"":" & Ref_channel_number(Current_Protocol) & ",")

        Print(1, """number_pulses"":[")
        Dim index As Integer
        For index = 1 To Number_Loops(Current_Protocol) - 1
            Print(1, Number_Pulses(Current_Protocol, index) & ",")

        Next
        Print(1, Number_Pulses(Current_Protocol, Number_Loops(Current_Protocol)) & "],")

        Print(1, """number_measuring_lights"":" & Number_Measuring_Lights(Current_Protocol) & ",")



        Print(1, """measuring_light"":[")
        'Dim index As Integer
        For index = 1 To Number_Measuring_Lights(Current_Protocol) - 1
            Print(1, Measuring_Light(Current_Protocol, index) & ",")

        Next
        Print(1, Measuring_Light(Current_Protocol, Number_Measuring_Lights(Current_Protocol)) & "],")

        'Measuring_Interval(, , )


        Print(1, """measuring_light_names"":[")
        'Dim index As Integer
        'Jabber_jabber.List1.Items.Add("measuring_light_names")
        For index = 0 To Number_Measuring_Lights(Current_Protocol) - 2
            Print(1, """" & M_Measuring_Light_Names(Current_Protocol, index) & """,")
            'Jabber_jabber.List1.Items.Add(M_Measuring_Light_Names(Current_Protocol, index))


        Next
        Print(1, """" & M_Measuring_Light_Names(Current_Protocol, Number_Measuring_Lights(Current_Protocol) - 1) & """],")
        Jabber_jabber.List1.Items.Add(M_Measuring_Light_Names(Current_Protocol, index))
        ' ,M_Measuring_Light_Names()

        'Intensity()
        Print(1, """actinic_intensity"":[")
        For index2 = 1 To Number_Loops(Current_Protocol)



            Print(1, Intensity(Current_Protocol, index2))

            If (index2 = Number_Loops(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next

        'M_Far_Red

        Print(1, """far_red"":[")
        For index2 = 1 To Number_Loops(Current_Protocol)



            Print(1, M_Far_Red(Current_Protocol, index2))

            If (index2 = Number_Loops(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next
        ' M_Blue_actinic

        Print(1, """blue_actinic"":[")
        For index2 = 1 To Number_Loops(Current_Protocol)



            Print(1, M_Blue_actinic(Current_Protocol, index2))

            If (index2 = Number_Loops(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next

        'Xe_Flash()

        Print(1, """flash"":[")
        For index2 = 1 To Number_Loops(Current_Protocol)



            Print(1, Xe_Flash(Current_Protocol, index2))

            If (index2 = Number_Loops(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next

        Print(1, """measuring_interval"":[")
        For index2 = 1 To Number_Loops(Current_Protocol)

            For index = 1 To Number_Measuring_Lights(Current_Protocol)
                Dim trimmed As String = Trim(Measuring_Interval(Current_Protocol, index, index2))
                If VB.Right(trimmed, 1) = "u" Then
                    factor = 0.000001
                ElseIf VB.Right(trimmed, 1) = "m" Then
                    factor = 0.001
                ElseIf VB.Right(trimmed, 1) = "n" Then
                    factor = 0.000000001
                Else
                    factor = 1
                End If
                Dim duration As Single = Val(trimmed) * factor

                Print(1, duration)

                If (index = Number_Measuring_Lights(Current_Protocol)) And (index2 = Number_Loops(Current_Protocol)) Then

                    Print(1, "],")
                Else

                    Print(1, ",")
                End If


            Next
        Next


        '        Measuring_Pulse_Duration(, )

        Print(1, """measuring_pulse_duration"":[")
        For index = 1 To Number_Measuring_Lights(Current_Protocol)
            Dim trimmed As String = Trim(Measuring_Pulse_Duration(Current_Protocol, index))
            If VB.Right(trimmed, 1) = "u" Then
                factor = 0.000001
            ElseIf VB.Right(trimmed, 1) = "m" Then
                factor = 0.001
            ElseIf VB.Right(trimmed, 1) = "n" Then
                factor = 0.000000001
            Else
                factor = 1
            End If
            Dim duration As Single = Val(trimmed) * factor

            Print(1, duration)

            If (index = Number_Measuring_Lights(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next


        '        Primary_Gain()
        Print(1, """sample_gain"":[")

        For index = 1 To Number_Measuring_Lights(Current_Protocol)
            Print(1, Primary_Gain(Current_Protocol, index))

            If (index = Number_Measuring_Lights(Current_Protocol)) Then

                Print(1, "],")
            Else

                Print(1, ",")
            End If


        Next
        'reference_gain

        Print(1, """reference_gain"":[")

        For index = 1 To Number_Measuring_Lights(Current_Protocol)
            Print(1, Primary_Gain(Current_Protocol, index))

            If (index = Number_Measuring_Lights(Current_Protocol)) Then

                Print(1, "]")
            Else

                Print(1, ",")
            End If


        Next
        Print(1, "}")
        FileClose(1)

        'Threading.Thread.sleep(10)

    End Sub
    Sub Plot_File(ByRef Plot_File_name As String, ByVal Plot_Window As Short, ByVal Plot_Delta As Short, ByVal add_to As Integer, ByVal linear_x As Integer, ByRef Current_Protocol As Short, ByVal m_plot_specific_traces(,) As Short) ', ByRef graph_name() As Object)

        Dim Number_Measuring_Lights As Integer
        Dim i, ii, total_points As Integer
        Dim j As String
        Dim Plot_Data_X(,) As Object
        Dim Plot_Data_Y(,) As Object
        Dim temp As Single = 0
        Dim strtemp As String = ""

        ' figure out how many columns
        FileOpen(10, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        s_Plot_File_Name(Plot_Window) = Plot_File_name
        s_Plot_File_Type(Plot_Window) = Plot_Delta

        j = LineInput(10)
        FileClose(10)

        '        List1.Items.Clear()
        'here we try to deduce the number of measuring lights
        'based on the fast that there will be n*4-1 tabs for each 
        ' measuring light (the n-1 is because I've eliminated the
        ' final (hanging) tab

        Number_Measuring_Lights = 1
        For i = 1 To Len(j)
            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))
            If Asc(Mid(j, i, 1)) = 9 Then  'counting the number of tabs
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        Dim Ttest As Integer = Number_Measuring_Lights

        Number_Measuring_Lights = Number_Measuring_Lights / 4
        FileOpen(11, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        ' in the following, we try to deduce the number of measuring lights
        ' based on the 
        i = 0
        'Dim strtemp As String
        While EOF(11) = False
            Input(11, temp)
            'temp = val(strtemp)
            'Input(11, temp)
            'Input(11, temp)
            'Input(11, temp)
            i = i + 1
        End While  'count the total number of points, each has 4 data
        FileClose(11)

        i = i / 4


        total_points = i / Number_Measuring_Lights
        ReDim Plot_Data_X(Number_Measuring_Lights, i)
        ReDim Plot_Data_Y(Number_Measuring_Lights, i)

        i = 0
        FileOpen(12, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        If Plot_Delta = 1 Then
            While EOF(12) = False

                For ii = 1 To Number_Measuring_Lights
                    Input(12, Plot_Data_X(ii, i))
                    If linear_x >= 1 Then
                        Plot_Data_X(ii, i) = i
                    End If
                    Input(12, temp)
                    Input(12, temp)
                    Input(12, Plot_Data_Y(ii, i))

                Next ii
                i = i + 1
            End While
        ElseIf Plot_Delta = 0 Then
            While EOF(12) = False
                For ii = 1 To Number_Measuring_Lights
                    Input(12, Plot_Data_X(ii, i))
                    Input(12, Plot_Data_Y(ii, i))
                    Input(12, temp)
                    Input(12, temp)

                Next ii
                i = i + 1
            End While
        ElseIf Plot_Delta = 2 Then  'plot reference
            While EOF(12) = False
                For ii = 1 To Number_Measuring_Lights
                    Input(12, Plot_Data_X(ii, i))
                    Input(12, temp)
                    Input(12, Plot_Data_Y(ii, i))
                    Input(12, temp)

                Next ii
                i = i + 1
            End While
        End If
        FileClose(12)


        Select Case Plot_Window
            Case 1
                Call CreateGraph(zg1, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(1), add_to, Plot_Window)
            Case 2
                Call CreateGraph(zg2, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(2), add_to, Plot_Window)
            Case 3
                Call CreateGraph(zg3, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(3), add_to, Plot_Window)
            Case 4
                Call CreateGraph(zg4, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(4), add_to, Plot_Window)
            Case 5
                Call CreateGraph(zg5, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(5), add_to, Plot_Window)
            Case 6
                Call CreateGraph(zg6, Number_Measuring_Lights, total_points, Plot_Data_X, Plot_Data_Y, Plot_Delta, Graph_Name(6), add_to, Plot_Window)

        End Select

    End Sub




    Private Sub UpdateLabels(ByRef This_Graph As Integer, ByRef Low_Val As Single, ByRef High_Val As Single)

    End Sub

    Sub Delta_A_calculate_smooth_reference(ByRef I_File_Name As String, ByVal baseline_start As Integer, ByVal baseline_end As Integer, ByRef smooth_window As Short, ByRef alternate_reference_file_name As String)
        Dim j As String
        Dim Number_Measuring_Lights As Integer

        'first, figure out how many measuring lights are in the file
        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        Dim strArray() As String
        strArray = j.Split(Chr(9))
        Jabber_jabber.List1.Items.Add(strArray.Length)
        Number_Measuring_Lights = strArray.Length / 4  '4 data points per measuring lights per 
        'next set up a matrix to hold the data sets
        Dim read_data_time(Number_Measuring_Lights, 1) As Single
        Dim read_data_sample(Number_Measuring_Lights, 1) As Single
        Dim read_data_reference(Number_Measuring_Lights, 1) As Single
        Dim read_data_reference_smoothed(Number_Measuring_Lights, 1) As Single

        Dim read_data_delta_a(Number_Measuring_Lights, 1) As Single


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        Dim i, ii As Integer
        ii = -1
        While EOF(1) = False
            ii = ii + 1
            ReDim Preserve read_data_time(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_sample(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_reference(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_reference_smoothed(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_delta_a(Number_Measuring_Lights, ii + 2)
            For i = 0 To Number_Measuring_Lights - 1
                Input(1, read_data_time(i, ii))
                Input(1, read_data_sample(i, ii))
                Input(1, read_data_reference(i, ii))
                Input(1, read_data_delta_a(i, ii))


            Next

        End While
        FileClose(1)
        Dim number_time_points As Integer = ii + 1

        If alternate_reference_file_name <> "" Then 'we are using a trace from anotehr run as an alternate reference 


            FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
            j = LineInput(1)
            FileClose(1)
            'Dim strArray() As String
            strArray = j.Split(Chr(9))
            Jabber_jabber.List1.Items.Add(strArray.Length)
            Number_Measuring_Lights = strArray.Length / 4  '4 data points per measuring lights per 
            'next set up a matrix to hold the data sets
            Dim a_read_data_time(Number_Measuring_Lights, 1) As Single
            Dim a_read_data_sample(Number_Measuring_Lights, 1) As Single
            Dim a_read_data_reference(Number_Measuring_Lights, 1) As Single
            Dim a_read_data_reference_smoothed(Number_Measuring_Lights, 1) As Single

            Dim a_read_data_delta_a(Number_Measuring_Lights, 1) As Single


            FileOpen(1, alternate_reference_file_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
            'Dim i, ii As Integer
            ii = -1
            While EOF(1) = False
                ii = ii + 1
                ReDim Preserve a_read_data_time(Number_Measuring_Lights, ii + 2)
                ReDim Preserve a_read_data_sample(Number_Measuring_Lights, ii + 2)
                ReDim Preserve a_read_data_reference(Number_Measuring_Lights, ii + 2)
                ReDim Preserve a_read_data_reference_smoothed(Number_Measuring_Lights, ii + 2)
                ReDim Preserve a_read_data_delta_a(Number_Measuring_Lights, ii + 2)
                For i = 0 To Number_Measuring_Lights - 1
                    Input(1, a_read_data_time(i, ii))
                    Input(1, a_read_data_sample(i, ii))
                    Input(1, a_read_data_reference(i, ii))
                    Input(1, a_read_data_delta_a(i, ii))


                Next

            End While
            FileClose(1)
            read_data_reference = a_read_data_sample 'use the sample data from one run as the refernece for another
        End If

        If smooth_window > 0 Then

            Dim ix, iix, count_avs As Integer
            Dim number_smooth_range As Integer = number_time_points

            For i = 0 To Number_Measuring_Lights - 1

                For ix = 0 To number_smooth_range - 1
                    count_avs = 0
                    If (ix >= smooth_window) And (ix < number_smooth_range - smooth_window) Then
                        For iix = ix - smooth_window To ix + smooth_window
                            read_data_reference_smoothed(i, ix) = read_data_reference_smoothed(i, ix) + read_data_reference(i, iix)
                            count_avs = count_avs + 1
                        Next
                    ElseIf (ix < smooth_window) Then
                        For iix = 0 To ix + smooth_window
                            read_data_reference_smoothed(i, ix) = read_data_reference_smoothed(i, ix) + read_data_reference(i, iix)
                            count_avs = count_avs + 1
                        Next
                    ElseIf ix >= number_smooth_range - smooth_window Then
                        For iix = ix - smooth_window To number_smooth_range - 1
                            read_data_reference_smoothed(i, ix) = read_data_reference_smoothed(i, ix) + read_data_reference(i, iix)

                            Jabber_jabber.List1.Items.Add(ix)
                            Jabber_jabber.List1.Items.Add(read_data_reference(i, ix))
                            count_avs = count_avs + 1
                        Next
                    End If
                    read_data_reference_smoothed(i, ix) = read_data_reference_smoothed(i, ix) / count_avs


                Next
            Next
        End If

        'calculate the baselines for each measuring light
        Dim baseline_offsets(Number_Measuring_Lights, 2) As Single
        Dim btemp1, btemp2 As Single
        btemp1 = btemp2 = 0
        If baseline_start > number_time_points Then
            baseline_start = 0
        End If

        If baseline_end > number_time_points Then
            baseline_end = number_time_points
        End If

        Dim sample_baselines(Number_Measuring_Lights) As Single
        Dim reference_baselines(Number_Measuring_Lights) As Single


        For i = 0 To Number_Measuring_Lights - 1
            btemp1 = 0
            btemp2 = 0
            For ii = baseline_start - 1 To baseline_end - 1
                btemp1 = btemp1 + read_data_sample(i, ii)
                If smooth_window > 0 Then

                    btemp2 = btemp2 + read_data_reference_smoothed(i, ii)
                Else
                    btemp2 = btemp2 + read_data_reference(i, ii)
                End If

            Next
            sample_baselines(i) = btemp1 / (baseline_end - baseline_start)
            reference_baselines(i) = btemp2 / (baseline_end - baseline_start)

        Next
        Dim xx(,) As Single = read_data_reference
        Dim x(,) As Single
        If smooth_window > 0 Then
            x = read_data_reference_smoothed
        Else
            x = read_data_reference
        End If
        'normalize traces to baseline averages
        Dim delta_a(Number_Measuring_Lights, number_time_points)
        Dim normalized_sample As Single
        Dim normalized_reference As Single

        For i = 0 To Number_Measuring_Lights - 1
            For ii = 0 To number_time_points - 1
                normalized_sample = read_data_sample(i, ii) / sample_baselines(i)
                If ref_mode = 1 Then
                    If smooth_window > 0 Then
                        normalized_reference = read_data_reference_smoothed(i, ii) / reference_baselines(i)
                    Else
                        normalized_reference = read_data_reference(i, ii) / reference_baselines(i)
                    End If
                Else
                    normalized_reference = 1
                End If
                delta_a(i, ii) = -1 * Math.Log(normalized_sample / normalized_reference)
            Next
        Next

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
        For ii = 0 To number_time_points - 1
            For i = 0 To Number_Measuring_Lights - 1
                Print(1, read_data_time(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_sample(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_reference(i, ii))
                Print(1, Chr(9))
                Print(1, delta_a(i, ii))
                If i < Number_Measuring_Lights - 1 Then
                    Print(1, Chr(9))
                End If
            Next
            PrintLine(1, "")
        Next
        FileClose(1)
    End Sub
    Sub Delta_A_calculate(ByRef I_File_Name As String, ByVal baseline_start As Integer, ByVal baseline_end As Integer)
        Dim j As String
        Dim Number_Measuring_Lights As Integer

        'first, figure out how many measuring lights are in the file
        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        Dim strArray() As String
        strArray = j.Split(Chr(9))
        'Jabber_jabber.List1.Items.Add(strArray.Length)
        Number_Measuring_Lights = strArray.Length / 4  '4 data points per measuring lights per 
        'next set up a matrix to hold the data sets
        Dim read_data_time(Number_Measuring_Lights, 1) As Single
        Dim read_data_sample(Number_Measuring_Lights, 1) As Single
        Dim read_data_reference(Number_Measuring_Lights, 1) As Single
        Dim read_data_delta_a(Number_Measuring_Lights, 1) As Single


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        Dim i, ii As Integer
        ii = -1
        While EOF(1) = False
            ii = ii + 1
            ReDim Preserve read_data_time(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_sample(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_reference(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_delta_a(Number_Measuring_Lights, ii + 2)
            For i = 0 To Number_Measuring_Lights - 1
                Input(1, read_data_time(i, ii))
                Input(1, read_data_sample(i, ii))
                Input(1, read_data_reference(i, ii))
                Input(1, read_data_delta_a(i, ii))


            Next

        End While
        FileClose(1)

        Dim number_time_points As Integer = ii + 1
        'calculate the baselines for each measuring light
        Dim baseline_offsets(Number_Measuring_Lights, 2) As Single
        Dim btemp1, btemp2 As Single
        btemp1 = btemp2 = 0
        If baseline_start > number_time_points Then
            baseline_start = 0
        End If

        If baseline_end > number_time_points Then
            baseline_end = number_time_points
        End If

        Dim sample_baselines(Number_Measuring_Lights) As Single
        Dim reference_baselines(Number_Measuring_Lights) As Single


        For i = 0 To Number_Measuring_Lights - 1
            btemp1 = 0
            btemp2 = 0
            For ii = baseline_start - 1 To baseline_end - 1
                btemp1 = btemp1 + read_data_sample(i, ii)
                btemp2 = btemp2 + read_data_reference(i, ii)
            Next
            sample_baselines(i) = btemp1 / (baseline_end - baseline_start)
            reference_baselines(i) = btemp2 / (baseline_end - baseline_start)

        Next

        'normalize traces to baseline averages
        Dim delta_a(Number_Measuring_Lights, number_time_points)
        Dim normalized_sample As Single
        Dim normalized_reference As Single

        For i = 0 To Number_Measuring_Lights - 1
            For ii = 0 To number_time_points - 1
                normalized_sample = read_data_sample(i, ii) / sample_baselines(i)
                If ref_mode = 1 Then
                    normalized_reference = read_data_reference(i, ii) / reference_baselines(i)
                Else
                    normalized_reference = 1
                End If
                delta_a(i, ii) = -1 * Math.Log(normalized_sample / normalized_reference)
            Next
        Next

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
        For ii = 0 To number_time_points - 1
            For i = 0 To Number_Measuring_Lights - 1
                Print(1, read_data_time(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_sample(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_reference(i, ii))
                Print(1, Chr(9))
                Print(1, delta_a(i, ii))
                If i < Number_Measuring_Lights - 1 Then
                    Print(1, Chr(9))
                End If
            Next
            PrintLine(1, "")
        Next
        FileClose(1)

    End Sub
    Sub normalize_to_integration_time(ByRef I_File_Name As String, ByVal Measuring_Pulse_Duration(,) As String, _
                                          ByRef Current_Protocol As Short)
        Dim j As String
        Dim Number_Measuring_Lights As Integer

        'first, figure out how many measuring lights are in the file
        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        Dim strArray() As String
        strArray = j.Split(Chr(9))
        Jabber_jabber.List1.Items.Add(strArray.Length)
        Number_Measuring_Lights = strArray.Length / 4  '4 data points per measuring lights per 
        'next set up a matrix to hold the data sets
        Dim read_data_time(Number_Measuring_Lights, 1) As Single
        Dim read_data_sample(Number_Measuring_Lights, 1) As Single
        Dim read_data_reference(Number_Measuring_Lights, 1) As Single
        Dim read_data_delta_a(Number_Measuring_Lights, 1) As Single


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        Dim i, ii As Integer
        ii = -1
        While EOF(1) = False
            ii = ii + 1
            ReDim Preserve read_data_time(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_sample(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_reference(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_delta_a(Number_Measuring_Lights, ii + 2)
            For i = 0 To Number_Measuring_Lights - 1
                Input(1, read_data_time(i, ii))
                Input(1, read_data_sample(i, ii))
                Input(1, read_data_reference(i, ii))
                Input(1, read_data_delta_a(i, ii))


            Next

        End While
        FileClose(1)

        Dim number_time_points As Integer = ii + 1
        'calculate the baselines for each measuring light
        'Dim baseline_offsets(Number_Measuring_Lights, 2) As Single
        'Dim btemp1, btemp2 As Single
        'btemp1 = btemp2 = 0
        'If baseline_start > number_time_points Then
        '    baseline_start = 0
        'End If

        'If baseline_end > number_time_points Then
        '    baseline_end = number_time_points
        'End If

        'Dim sample_baselines(Number_Measuring_Lights) As Single
        'Dim reference_baselines(Number_Measuring_Lights) As Single


        'For i = 0 To Number_Measuring_Lights - 1
        '    btemp1 = 0
        '    btemp2 = 0
        '    For ii = baseline_start - 1 To baseline_end - 1
        '        btemp1 = btemp1 + read_data_sample(i, ii)
        '        btemp2 = btemp2 + read_data_reference(i, ii)
        '    Next
        '    sample_baselines(i) = btemp1 / (baseline_end - baseline_start)
        '    reference_baselines(i) = btemp2 / (baseline_end - baseline_start)

        'Next

        'normalize traces to baseline averages
        Dim normalized_values(Number_Measuring_Lights, number_time_points)
        Dim normalized_sample As Single
        Dim normalized_reference As Single

        Dim integration_time As Single
        Jabber_jabber.List1.Items.Add("normalize to integration time")
        Dim trimmed As String
        Dim factor As Double
        Dim duration As Double
        For i = 0 To Number_Measuring_Lights - 1
            trimmed = Trim(Measuring_Pulse_Duration(Current_Protocol, i + 1))
            If VB.Right(trimmed, 1) = "u" Then
                factor = 0.000001
            ElseIf VB.Right(trimmed, 1) = "m" Then
                factor = 0.001
            ElseIf VB.Right(trimmed, 1) = "n" Then
                factor = 0.000000001
            Else
                factor = 1
            End If
            integration_time = Val(trimmed) * factor

            Jabber_jabber.List1.Items.Add(integration_time)
            For ii = 0 To number_time_points - 1
                normalized_sample = read_data_sample(i, ii) / integration_time
                Jabber_jabber.List1.Items.Add(normalized_sample)
                'If ref_mode = 1 Then
                normalized_reference = read_data_reference(i, ii) / integration_time
                'Else
                'normalized_reference = 1
                'End If
                normalized_values(i, ii) = normalized_sample
            Next
        Next

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
        For ii = 0 To number_time_points - 1
            For i = 0 To Number_Measuring_Lights - 1
                Print(1, read_data_time(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_sample(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_reference(i, ii))
                Print(1, Chr(9))
                Print(1, normalized_values(i, ii))
                If i < Number_Measuring_Lights - 1 Then
                    Print(1, Chr(9))
                End If
            Next
            PrintLine(1, "")
        Next
        FileClose(1)

    End Sub

    Sub Ratio_Calculate(ByRef I_File_Name As String, ByVal baseline_start As Integer, ByVal baseline_end As Integer)
        Dim j As String
        Dim Number_Measuring_Lights As Integer

        'first, figure out how many measuring lights are in the file
        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        Dim strArray() As String
        strArray = j.Split(Chr(9))
        Jabber_jabber.List1.Items.Add(strArray.Length)
        Number_Measuring_Lights = strArray.Length / 4  '4 data points per measuring lights per 
        'next set up a matrix to hold the data sets
        Dim read_data_time(Number_Measuring_Lights, 1) As Single
        Dim read_data_sample(Number_Measuring_Lights, 1) As Single
        Dim read_data_reference(Number_Measuring_Lights, 1) As Single
        Dim read_data_delta_a(Number_Measuring_Lights, 1) As Single


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        Dim i, ii As Integer
        ii = -1
        While EOF(1) = False
            ii = ii + 1
            ReDim Preserve read_data_time(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_sample(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_reference(Number_Measuring_Lights, ii + 2)
            ReDim Preserve read_data_delta_a(Number_Measuring_Lights, ii + 2)
            For i = 0 To Number_Measuring_Lights - 1
                Input(1, read_data_time(i, ii))
                Input(1, read_data_sample(i, ii))
                Input(1, read_data_reference(i, ii))
                Input(1, read_data_delta_a(i, ii))


            Next

        End While
        FileClose(1)

        Dim number_time_points As Integer = ii + 1
        'calculate the baselines for each measuring light
        Dim baseline_offsets(Number_Measuring_Lights, 2) As Single
        Dim btemp1, btemp2 As Single
        btemp1 = btemp2 = 0
        If baseline_start > number_time_points Then
            baseline_start = 0
        End If

        If baseline_end > number_time_points Then
            baseline_end = number_time_points
        End If

        Dim sample_baselines(Number_Measuring_Lights) As Single
        Dim reference_baselines(Number_Measuring_Lights) As Single


        For i = 0 To Number_Measuring_Lights - 1
            btemp1 = 0
            btemp2 = 0
            For ii = baseline_start - 1 To baseline_end - 1
                btemp1 = btemp1 + read_data_sample(i, ii)
                btemp2 = btemp2 + read_data_reference(i, ii)
            Next
            sample_baselines(i) = btemp1 / (baseline_end - baseline_start)
            reference_baselines(i) = btemp2 / (baseline_end - baseline_start)

        Next

        'normalize traces to baseline averages
        Dim ratio_values(Number_Measuring_Lights, number_time_points)
        Dim normalized_sample As Single
        Dim normalized_reference As Single

        For i = 0 To Number_Measuring_Lights - 1
            For ii = 0 To number_time_points - 1
                normalized_sample = read_data_sample(i, ii) / sample_baselines(i)
                'If ref_mode = 1 Then
                normalized_reference = read_data_reference(i, ii) / reference_baselines(i)
                'Else
                'normalized_reference = 1
                'End If
                ratio_values(i, ii) = -1 * Math.Log(normalized_sample / normalized_reference)
            Next
        Next

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
        For ii = 0 To number_time_points - 1
            For i = 0 To Number_Measuring_Lights - 1
                Print(1, read_data_time(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_sample(i, ii))
                Print(1, Chr(9))
                Print(1, read_data_reference(i, ii))
                Print(1, Chr(9))
                Print(1, ratio_values(i, ii))
                If i < Number_Measuring_Lights - 1 Then
                    Print(1, Chr(9))
                End If
            Next
            PrintLine(1, "")
        Next
        FileClose(1)

    End Sub
    'Sub Delta_A_no_ref_old(ByRef I_File_Name As String, ByVal baseline_start As Integer, ByVal baseline_end As Integer)

    '    Dim Number_columns As Object
    '    Dim Number_Measuring_Lights As Short
    '    Dim i, ii As Short
    '    Dim total_points As Short
    '    Dim j As String
    '    Dim p, pp As Integer
    '    'Dim ppp, p, pp, light As Object
    '    'Dim Column_N As Short

    '    Dim Delta_Data_X() As Object
    '    Dim Delta_Data_Y() As Object
    '    Dim Delta_Data_ref() As Object
    '    Dim temp As Single = 0
    '    Dim tempp As Single
    '    'Dim I0_Y As Single
    '    Dim I0_Data_Y() As Single
    '    Dim av_baseline() As Single

    '    ' figure out how many columns are in I_FILE


    '    FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
    '    j = LineInput(1)
    '    FileClose(1)
    '    '        Threading.Thread.sleep(10)



    '    Number_Measuring_Lights = 1
    '    For i = 1 To Len(j)

    '        'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

    '        If Asc(Mid(j, i, 1)) = 9 Then
    '            Number_Measuring_Lights = Number_Measuring_Lights + 1
    '        End If
    '    Next i

    '    'List1.AddItem j$
    '    'List1.AddItem "l: " & Len(j$)

    '    Number_Measuring_Lights = Number_Measuring_Lights / 4

    '    Number_columns = Number_Measuring_Lights


    '    FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

    '    i = 0

    '    While EOF(1) = False
    '        Input(1, temp)
    '        Input(1, temp)
    '        Input(1, temp)
    '        Input(1, temp)
    '        i = i + 1
    '    End While
    '    FileClose(1)
    '    'Threading.Thread.sleep(10)

    '    total_points = i
    '    'dimension the temporary arrays

    '    ReDim Delta_Data_X(i)
    '    ReDim Delta_Data_Y(i)
    '    ReDim Delta_Data_ref(i)

    '    ReDim I0_Data_Y(i)
    '    ReDim av_baseline(Number_Measuring_Lights)

    '    ' input the RAW data from I_FILE

    '    i = 0

    '    FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)


    '    While EOF(1) = False
    '        Input(1, Delta_Data_X(i))
    '        Input(1, Delta_Data_Y(i))
    '        Input(1, Delta_Data_ref(i))
    '        Input(1, temp)
    '        i = i + 1
    '    End While

    '    FileClose(1)
    '    'Threading.Thread.sleep(10)

    '    For pp = 0 To Number_Measuring_Lights - 1
    '        av_baseline(pp) = 0
    '    Next pp
    '    i = 0
    '    ii = 0
    '    For p = 0 To total_points - Number_Measuring_Lights Step Number_Measuring_Lights

    '        If i >= baseline_start And i <= baseline_end Then
    '            ii = ii + 1
    '            For pp = 0 To Number_Measuring_Lights - 1
    '                av_baseline(pp) = av_baseline(pp) + Delta_Data_Y(i + pp)
    '            Next pp
    '        End If
    '        i = i + 1
    '    Next p

    '    For pp = 0 To Number_Measuring_Lights - 1
    '        av_baseline(pp) = av_baseline(pp) / ii
    '        List1.Items.Add("baseline=" + Str(pp + 1) + "-->" + Str(av_baseline(pp)))

    '    Next pp
    '    ' input the RAW data from I0_FILE

    '    i = 0


    '    ' now open I0_FILE, calculate (deltaA, then resave the data into I_FILE

    '    FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

    '    For i = 0 To total_points - Number_columns Step Number_columns
    '        'Column_N = 0


    '        For p = 0 To Number_Measuring_Lights - 1

    '            If (Delta_Data_Y(i + p) / av_baseline(p)) > 0 Then
    '                tempp = -1 * (Math.Log10(Delta_Data_Y(i + p) / av_baseline(p)))
    '            Else
    '                tempp = 0
    '            End If
    '            If p < Number_Measuring_Lights Then
    '                Print(1, Delta_Data_X(i + p) & _
    '                                      Chr(9) & Delta_Data_Y(i + p) & _
    '                                      Chr(9) & Delta_Data_ref(i + p) & _
    '                                      Chr(9) & tempp & Chr(9))
    '            Else
    '                Print(1, Delta_Data_X(i + p) & _
    '                                      Chr(9) & Delta_Data_Y(i + p) & _
    '                                      Chr(9) & Delta_Data_ref(i + p) & _
    '                                      Chr(9) & tempp)
    '            End If
    '        Next p
    '        PrintLine(1, "")
    '    Next i
    '    FileClose(1)
    '    'Threading.Thread.sleep(10)


    'End Sub

    Sub Delta_A_WITH_ref(ByRef I_File_Name As String, ByVal baseline_start As Integer, ByVal baseline_end As Integer)

        Dim Number_columns As Object
        Dim Number_Measuring_Lights As Short
        Dim i, ii As Short
        Dim total_points As Short
        Dim j As String
        Dim p, pp As Integer
        'Dim ppp, p, pp, light As Object
        'Dim Column_N As Short

        Dim Delta_Data_X() As Object
        Dim Delta_Data_Y() As Object
        Dim Delta_Data_ref() As Object
        Dim temp As Single = 0
        Dim tempp As Single
        'Dim I0_Y As Single
        Dim I0_Data_Y() As Single

        Dim av_baseline() As Single

        ' figure out how many columns are in I_FILE


        FileOpen(13, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(13)
        FileClose(13)
        'Threading.Thread.sleep(10)




        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        'List1.AddItem j$
        'List1.AddItem "l: " & Len(j$)

        Number_Measuring_Lights = Number_Measuring_Lights / 4

        Number_columns = Number_Measuring_Lights


        FileOpen(14, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0

        While EOF(14) = False
            Input(14, temp)
            Input(14, temp)
            Input(14, temp)
            Input(14, temp)
            i = i + 1
        End While
        FileClose(14)
        'Threading.Thread.sleep(10)

        total_points = i
        'dimension the temporary arrays

        ReDim Delta_Data_X(i)
        ReDim Delta_Data_Y(i)
        ReDim Delta_Data_ref(i)

        ReDim I0_Data_Y(i)
        ReDim av_baseline(Number_Measuring_Lights)

        ' input the RAW data from I_FILE

        i = 0

        FileOpen(15, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)


        While EOF(15) = False
            Input(15, Delta_Data_X(i))
            Input(15, Delta_Data_Y(i))
            Input(15, Delta_Data_ref(i))
            Input(15, temp)
            i = i + 1
        End While

        FileClose(15)
        'Threading.Thread.sleep(10)

        For pp = 0 To Number_Measuring_Lights - 1
            av_baseline(pp) = 0
        Next pp
        i = 0
        ii = 0
        For p = 0 To total_points - Number_Measuring_Lights Step Number_Measuring_Lights

            If i >= baseline_start And i <= baseline_end Then
                ii = ii + 1
                For pp = 0 To Number_Measuring_Lights - 1

                    If Delta_Data_ref(p + pp) <> 0 Then

                        If (Delta_Data_Y(p + pp) / Delta_Data_ref(p + pp)) > 0 Then
                            av_baseline(pp) = av_baseline(pp) + (Math.Log10(Delta_Data_ref(p + pp) / Delta_Data_Y(p + pp)))
                        End If
                    End If

                Next pp
            End If
            i = i + 1
        Next p

        For pp = 0 To Number_Measuring_Lights - 1
            av_baseline(pp) = av_baseline(pp) / ii
        Next pp
        ' input the RAW data from I0_FILE

        i = 0


        ' now open I0_FILE, calculate (deltaA, then resave the data into I_FILE

        FileOpen(16, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        For i = 0 To total_points - Number_columns Step Number_columns
            'Column_N = 0


            For p = 0 To Number_Measuring_Lights - 1
                If Delta_Data_ref(i + p) <> 0 Then
                    If (Delta_Data_Y(i + p) / Delta_Data_ref(i + p)) > 0 Then
                        tempp = (Math.Log10(Delta_Data_ref(i + p) / Delta_Data_Y(i + p))) - av_baseline(p)
                    Else
                        tempp = 0
                    End If
                End If
                If p < Number_Measuring_Lights Then
                    Print(1, Delta_Data_X(i + p) & _
                                          Chr(9) & Delta_Data_Y(i + p) & _
                                          Chr(9) & Delta_Data_ref(i + p) & _
                                          Chr(9) & tempp & Chr(9))
                Else
                    Print(1, Delta_Data_X(i + p) & _
                                          Chr(9) & Delta_Data_Y(i + p) & _
                                          Chr(9) & Delta_Data_ref(i + p) & _
                                          Chr(9) & tempp)
                End If

            Next p
            PrintLine(16, "")
        Next i
        FileClose(16)

        'Threading.Thread.sleep(10)

    End Sub

    Sub Zero_Point(ByRef I_File_Name As String, ByVal point_to_zero As Integer)

        Dim Number_columns As Object
        Dim Number_Measuring_Lights As Short
        Dim i, ii As Short
        Dim total_points As Short
        Dim j As String
        Dim p, pp As Integer
        'Dim ppp, p, pp, light As Object
        'Dim Column_N As Short

        Dim Delta_Data_X() As Object
        Dim Delta_Data_Y() As Object
        Dim Delta_Data_ref() As Object
        Dim Delta_Data_Delta() As Object
        Dim temp As Single = 0
        Dim tempp As Single
        'Dim I0_Y As Single
        Dim I0_Data_Y() As Single
        Dim av_baseline() As Single

        ' figure out how many columns are in I_FILE


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        'List1.AddItem j$
        'List1.AddItem "l: " & Len(j$)

        Number_Measuring_Lights = Number_Measuring_Lights / 4

        Number_columns = Number_Measuring_Lights


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0

        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        total_points = i

        'dimension the temporary arrays

        ReDim Delta_Data_X(i)
        ReDim Delta_Data_Y(i)
        ReDim Delta_Data_ref(i)
        ReDim Delta_Data_Delta(i)
        ReDim I0_Data_Y(i)
        ReDim av_baseline(Number_Measuring_Lights)

        ' input the RAW data from I_FILE

        i = 0

        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)


        While EOF(1) = False
            Input(1, Delta_Data_X(i))
            Input(1, Delta_Data_Y(i))
            Input(1, Delta_Data_ref(i))
            Input(1, Delta_Data_Delta(i))
            i = i + 1
        End While

        FileClose(1)
        'Threading.Thread.sleep(10)


        i = 0
        ii = 0
        For p = 0 To total_points - Number_Measuring_Lights Step Number_Measuring_Lights
            i = i + 1
            If i = point_to_zero Then

                For pp = 0 To Number_Measuring_Lights - 1


                Next pp
            End If
            i = i + 1
        Next p

        ' input the RAW data from I0_FILE

        i = 0


        ' now open I0_FILE, calculate (deltaA, then resave the data into I_FILE

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        For i = 0 To total_points - Number_columns Step Number_columns
            'Column_N = 0


            For p = 0 To Number_Measuring_Lights - 1

                If i = point_to_zero - 1 Then
                    tempp = 0
                Else
                    tempp = Delta_Data_Delta(i + p)
                End If

                Print(1, Delta_Data_X(i + p) & Chr(9) & Delta_Data_Y(i + p) & Chr(9) & Delta_Data_ref(i + p) & Chr(9) & tempp & Chr(9))
            Next p
            PrintLine(1, "")
        Next i
        FileClose(1)

        ' Threading.Thread.sleep(10)

    End Sub



    Sub delete_point(ByRef I_File_Name As String, ByVal point_to_delete As Integer)

        Dim Number_columns As Object
        Dim Number_Measuring_Lights As Short
        Dim i As Short
        'Dim i, ii As Short
        Dim total_points As Short
        Dim j As String
        'Dim p, pp As Integer
        Dim p As Integer
        'Dim ppp, p, pp, light As Object
        'Dim Column_N As Short

        Dim Delta_Data_X() As Object
        Dim Delta_Data_Y() As Object
        Dim Delta_Data_ref() As Object
        Dim temp As Single = 0
        Dim tempp As Single
        'Dim I0_Y As Single
        Dim I0_Data_Y() As Single
        Dim av_baseline() As Single

        ' figure out how many columns are in I_FILE


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        'List1.AddItem j$
        'List1.AddItem "l: " & Len(j$)

        Number_Measuring_Lights = Number_Measuring_Lights / 4

        Number_columns = Number_Measuring_Lights


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0

        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        total_points = i
        'dimension the temporary arrays

        ReDim Delta_Data_X(i)
        ReDim Delta_Data_Y(i)
        ReDim Delta_Data_ref(i)

        ReDim I0_Data_Y(i)
        ReDim av_baseline(Number_Measuring_Lights)

        ' input the RAW data from I_FILE

        i = 0

        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)


        While EOF(1) = False
            Input(1, Delta_Data_X(i))
            Input(1, Delta_Data_Y(i))
            Input(1, Delta_Data_ref(i))
            Input(1, temp)
            i = i + 1
        End While

        FileClose(1)

        'Threading.Thread.sleep(10)


        ' now open I0_FILE, then resave the data into I_FILE

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        For i = 0 To total_points - Number_columns Step Number_columns
            'Column_N = 0


            For p = 0 To Number_Measuring_Lights - 1

                If (Delta_Data_Y(i + p) / av_baseline(p)) > 0 Then
                    tempp = (Math.Log10(Delta_Data_Y(i + p) / av_baseline(p)))
                Else
                    tempp = 0
                End If
                Print(1, Delta_Data_X(i + p) & Chr(9) & Delta_Data_Y(i + p) & Chr(9) & Delta_Data_ref(i + p) & Chr(9) & tempp)
            Next p
            PrintLine(1, "")
        Next i
        FileClose(1)
        'Threading.Thread.sleep(10)


    End Sub

    Sub Delta(ByRef I_File_Name As String, ByRef I0_File_Name As String)

        Dim Number_columns As Object
        Dim Number_Measuring_Lights As Short
        Dim i As Short
        Dim j As String
        Dim p, pp As Integer
        'Dim ppp, p, pp, light As Object
        'Dim Column_N As Short

        Dim Delta_Data_X() As Object
        Dim Delta_Data_Y() As Object
        Dim temp As Single = 0
        'Dim I0_Y As Single
        Dim I0_Data_Y() As Single

        ' figure out how many columns are in I_FILE


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        ' Threading.Thread.sleep(10)

        '        List1.Items.Clear()

        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        'List1.AddItem j$
        'List1.AddItem "l: " & Len(j$)

        Number_Measuring_Lights = Number_Measuring_Lights / 4

        Number_columns = Number_Measuring_Lights


        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0

        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        'dimension the temporary arrays

        ReDim Delta_Data_X(i)
        ReDim Delta_Data_Y(i)
        ReDim I0_Data_Y(i)

        ' input the RAW data from I_FILE

        i = 0

        FileOpen(1, I_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        While EOF(1) = False
            Input(1, Delta_Data_X(i))
            Input(1, Delta_Data_Y(i))
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While

        FileClose(1)

        'Threading.Thread.sleep(10)

        ' input the RAW data from I0_FILE

        i = 0

        FileOpen(1, I0_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        While EOF(1) = False
            Input(1, temp)
            Input(1, I0_Data_Y(i))
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While

        FileClose(1)
        'Threading.Thread.sleep(10)

        ' now open I0_FILE, calculate (I-I0)/I0, then resave the data into I_FILE

        FileOpen(1, I_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)



        For pp = 0 To i - 1 Step Number_columns
            'Column_N = 0


            For p = pp To pp + Number_columns - 1

                temp = Int(1000000 * ((Delta_Data_Y(p) / I0_Data_Y(p)) - 1)) / 1000000

                If p > pp Then Print(1, Chr(9))


                Print(1, Delta_Data_X(p) & Chr(9) & Int(1000000 * Delta_Data_Y(p)) / 1000000 & Chr(9) & temp & Chr(9) & temp)
            Next p
            PrintLine(1, "")
        Next pp
        FileClose(1)

        'Threading.Thread.sleep(10)

    End Sub



    Sub Dirk_Save(ByRef data_time() As Single, ByRef points As Integer, ByRef File_Name As String, _
                  ByRef Save_Mode As Short, ByRef Times_Averaged As Short, ByRef S_Time_Mode As String, _
                  ByRef Start_Time As Single, ByRef Number_columns As Short, ByRef Baseline_Start As Object, _
                  ByRef Baseline_End As Object, ByVal In_Channel_in As Short, ByVal Ref_Channel_number_in As Short, _
                  ByRef invert_raw As Integer, ByRef take_offset As Boolean)

        Jabber_jabber.FilesListBox.Items.Add(File_Name)
        Dim p, pp As Object
        Dim ii As Integer
        Dim Column_N As Short
        Dim Time_Offset As Single
        ' Dim ref_channel = 2
        'Dim j As String
        Dim junk As Single
        'Dim Reference_Data_Y() As Single

        Dim Accumulated_Data_X() As Single
        Dim Accumulated_Data_Y() As Single
        Dim Accumulated_Data_Ref() As Single

        ReDim Accumulated_Data_X(points)
        ReDim Accumulated_Data_Y(points)
        ReDim Accumulated_Data_Ref(points)
        Dim Reference_Data_Y(points) As Single

        'ReDim Baseline_Y(Number_columns)
        'remed out put back:



        If record_files = True Then


            note_text = note_text & vbCrLf & "'tod':" & "'" & TimeOfDay() & "'," & "'file':" & "'" & File_Name & "'" & vbCrLf

            FileOpen(7, Note_File_Name, OpenMode.Output)
            Write(7, note_text)
            FileClose(7)
        End If
        '...
        System.Windows.Forms.Application.DoEvents()

        If Baseline_Start > points Then Baseline_Start = 1 'prevents overflow errors when baseline points exceed the range

        If Baseline_End > points Then Baseline_End = points

        's
        ' take offset from all raw values
        ' NOTE: Dava_Volts is 0-indexed

        If take_offset = 1 Then
            For pp = 0 To points - 1 'data_volts is zero indexed
                For ii = 0 To 3
                    If ii = In_Channel_in Or ii = Ref_Channel_number_in Then
                        Data_Volts(pp, ii) = invert_raw * (Data_Volts(pp, ii) - offset)
                    End If

                Next ii
            Next pp
        End If

        If File_Name = "" Then
            Jabber_jabber.FilesListBox.Items.Add("no file name!")
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1
            Exit Sub
        End If

        If S_Time_Mode = "from_zero" Then
            Time_Offset = 0
        ElseIf S_Time_Mode = "sequential" Then
            Time_Offset = Start_Time
        End If

        '        Save_Mode = File_Replace

        Dim t As Integer = Data_Volts.Length()

        If Save_Mode = File_Replace Then
        End If 

        If Save_Mode = File_Append Then

            ' calcluate baselines
            If File.Exists(File_Name) = False Then 'if the file already exists use previous values for baseline
                For pp = Baseline_Start * Number_columns To Baseline_End * Number_columns Step Number_columns
                    Column_N = 0
                    For p = pp To pp + Number_columns - 1
                        Column_N = Column_N + 1
                        Baseline_Y(Column_N) = Baseline_Y(Column_N) + Data_Volts(p, In_Channel_in)
                    Next p
                Next pp


                For Column_N = 1 To Number_columns

                    Baseline_Y(Column_N) = Baseline_Y(Column_N) / (Baseline_End - Baseline_Start + 1)

                Next Column_N
            End If

            FileClose(1)
            FileOpen(1, File_Name, OpenMode.Append, OpenAccess.Default, OpenShare.Shared)

            For pp = 0 To points - 1 Step Number_columns 'here, we are assuming that data_volts is zero indexed
                '                                           and the Number_columns is equal to the number of lights???
                Column_N = 0

                For p = pp To pp + Number_columns - 1
                    Column_N = Column_N + 1
                    If p > pp Then Print(1, Chr(9))

                    '***************************WITH EXTRA CHANNEL************************

                    Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel_in) & Chr(9))
                    'If extra_channel = True Then
                    ' Print(1, (Data_Volts(p, 1)) & Chr(9)) '*************************WITH EXTRA CHANNEL*********************************
                    ' End If

                    Print(1, (Data_Volts(p, Ref_Channel_number_in)) & Chr(9))

                    If ref_mode = 1 Then

                        If Math.Abs(Data_Volts(p, Ref_Channel_number_in)) > 0 And (Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)) > 0 Then
                            Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)))
                        Else
                            Print(1, 0)
                        End If
                    Else

                        'Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel) & Chr(9))
                        'List1.AddItem p & " " & data_volts(p)
                        '                If Baseline_Y(Column_N) <> 0 Then
                        'Print(1, (Data_Volts(p, In_Channel) / Data_Volts(p, In_Channel + 1) - 1))
                        'Print(1, ((Data_Volts(p, In_Channel) / Baseline_Y(Column_N)) - 1))
                        'Print #1, (Int(1000000 * (Data_Volts(p) / Baseline_Y(Column_N)) - 1)) / 1000000;
                        'Else
                        'Print(1, 0)
                        'End If
                        'List1.AddItem p & " " & data_volts(p)
                        If Math.Abs(Baseline_Y(Column_N)) > 0 And (Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N)) > 0 Then
                            Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N)))
                        Else
                            Print(1, 0)
                        End If
                    End If
                Next p
                PrintLine(1, "") 'equivalent to CRLF??
            Next pp
            FileClose(1)




        ElseIf (Save_Mode = Average_Into) or (Save_Mode = Sequential) Then
            If Times_Averaged = 1 Then ' first time off, just save trace

                ' calcluate baselines


                For pp = Baseline_Start * Number_columns To Baseline_End * Number_columns Step Number_columns
                    Column_N = 0
                    For p = pp To pp + Number_columns - 1
                        Column_N = Column_N + 1
                        Baseline_Y(Column_N) = Baseline_Y(Column_N) + Data_Volts(p, In_Channel_in)
                    Next p
                Next pp

                For Column_N = 1 To Number_columns
                    Baseline_Y(Column_N) = Baseline_Y(Column_N) / (Baseline_End - Baseline_Start + 1)

                    'List1.AddItem "col " & Column_N & " " & Baseline_Y(Column_N)
                Next Column_N

                FileOpen(1, File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
                '       Print #1, "m", Number_Columns
                For pp = 0 To points - 1 Step Number_columns
                    Column_N = 0

                    For p = pp To pp + Number_columns - 1
                        Column_N = Column_N + 1
                        'If Column_N > 1 Then Print(1, Chr(9))
                        If p > pp Then Print(1, Chr(9))

                        Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel_in) & Chr(9))
                        Print(1, (Data_Volts(p, Ref_Channel_number_in)) & Chr(9))

                        If ref_mode = 1 Then
                            If Math.Abs(Data_Volts(p, Ref_Channel_number_in)) > 0 And (Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)) > 0 Then
                                Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)))
                            Else
                                Print(1, 0)
                            End If
                        Else
                            If Math.Abs(Baseline_Y(Column_N)) > 0 And (Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N)) > 0 Then
                                Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N)))
                            Else
                                Print(1, 0)
                            End If

                        End If


                        '                        Print(1, ((Data_Volts(p, In_Channel) / Baseline_Y(Column_N)) - 1))

                        '                        Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel) & Chr(9))

                        '                       Print(1, (Data_Volts(p, In_Channel) / Data_Volts(p, In_Channel + 1) - 1))

                        'Print(1, Int(1000000 * ((Data_Volts(p) / Baseline_Y(Column_N)) - 1)) / 1000000)
                        'Print #1, (Int(1000000 * (Data_Volts(p) / Baseline_Y(Column_N)) - 1)) / 1000000;

                        'Else
                        'Print(1, 0)
                        'End If
                        'List1.AddItem p & " " & data_volts(p)
                    Next p
                    PrintLine(1, "")
                Next pp
                FileClose(1)
                'Threading.Thread.sleep(10)


            Else
                FileOpen(1, File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data
                For pp = 0 To points - 1 Step Number_columns

                    For p = pp To pp + Number_columns - 1
                        Input(1, junk)
                        Input(1, Accumulated_Data_Y(p))
                        Input(1, Accumulated_Data_Ref(p))
                        Input(1, junk) ' throw away the deltaI/I data
                    Next p

                Next pp
                FileClose(1)
                'Threading.Thread.sleep(10)

                'next make average with current data

                For p = 0 To points - 1
                    Data_Volts(p, In_Channel_in) = (Data_Volts(p, In_Channel_in) + ((Times_Averaged - 1) * Accumulated_Data_Y(p))) / Times_Averaged
                    Data_Volts(p, Ref_Channel_number_in) = (Data_Volts(p, Ref_Channel_number_in) + ((Times_Averaged - 1) * Accumulated_Data_Ref(p))) / Times_Averaged
                Next p


                ' calcluate baselines

                For pp = Baseline_Start * Number_columns To Baseline_End * Number_columns Step Number_columns
                    Column_N = 0
                    For p = pp To pp + Number_columns - 1
                        Column_N = Column_N + 1
                        Baseline_Y(Column_N) = Baseline_Y(Column_N) + Data_Volts(p, In_Channel_in)
                    Next p
                Next pp

                For Column_N = 1 To Number_columns
                    Baseline_Y(Column_N) = Baseline_Y(Column_N) / (Baseline_End - Baseline_Start + 1)
                Next Column_N



                FileOpen(1, File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
                For pp = 0 To points - 1 Step Number_columns
                    Column_N = 0
                    For p = pp To pp + Number_columns - 1
                        Column_N = Column_N + 1
                        If p > pp Then Print(1, Chr(9))


                        Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel_in) & Chr(9))
                        Print(1, (Data_Volts(p, Ref_Channel_number_in)) & Chr(9))
                        'Print(1, Math.Log10(Data_Volts(p, In_Channel) / Data_Volts(p, 2)))

                        If ref_mode = 1 Then


                            If Math.Abs(Data_Volts(p, Ref_Channel_number_in)) > 0 And (Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)) > 0 Then
                                Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)))
                            Else
                                Print(1, 0)
                            End If
                        Else  'ref_channel is zero

                            If Math.Abs(Data_Volts(p, Ref_Channel_number_in)) > 0 And (Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N) > 0) And (Data_Volts(p, In_Channel_in) / Data_Volts(p, Ref_Channel_number_in)) > 0 Then
                                Print(1, -1 * Math.Log10(Data_Volts(p, In_Channel_in) / Baseline_Y(Column_N)))
                            Else
                                Print(1, 0)
                            End If
                        End If
                        '                        Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel) & Chr(9))
                        'List1.AddItem p & " " & data_volts(p)
                        'If Baseline_Y(Column_N) <> 0 Then
                        'Print #1, (Data_Volts(p) / Baseline_Y(Column_N)) - 1;
                        '                       Print(1, (Data_Volts(p, In_Channel) / Data_Volts(p, In_Channel + 1) - 1))
                        'Print(1, Int(1000000 * ((Data_Volts(p) / Baseline_Y(Column_N)) - 1)) / 1000000)
                        'Print #1, (Int(1000000 * (Data_Volts(p) / Baseline_Y(Column_N)) - 1)) / 1000000;
                        'Else
                        'Print(1, 0)
                        'End If
                        'List1.AddItem p & " " & data_volts(p)
                    Next p
                    PrintLine(1, "")
                Next pp
                FileClose(1)

            End If
        End If


    End Sub

    Function Auto_Gain_Calc(ByRef Current_Protocol As Integer, ByRef data_time() As Single, ByRef points As Integer, ByRef File_Name As String, ByRef Number_columns As Short, ByVal _
                            In_Channel() As Short, ByRef Ref_Channel_Number() As Short, ByRef length_of_duration_sequence As Integer, _
                            ByVal Measuring_Pulse_Duration(,) As String, _
                            ByRef exclude_from_gain As Integer, ByRef invert_raw As Single, _
                            ByRef primary_gain_index As Integer)
        'Dim voltages = {}


        Dim pp As Object
        Dim i, ii As Integer
        'Dim Column_N As Short
        Dim Upper_Value As Integer = 0

        'Dim Gain_Volts(2, 32) As Single
        'Dim m, b, dur1, dur2 As Single
        'Dim target As Single
        'Dim target_string As String



        For pp = 0 To points - 1
            For ii = 0 To 3
                Data_Volts(pp, ii) = invert_raw * (Data_Volts(pp, ii) - offset)
            Next ii
        Next pp
        Dim indexit As Integer
        i = 0
        'gain_data_list.Items.Add("-----------primary_gain = " & Str(i))
        For ii = 0 To length_of_duration_sequence - 1
            indexit = ii + (i * (length_of_duration_sequence))
            'gain_data_list.Items.Add(Str(indexit) & " pg:" & Str(i) & " t:" & Str(data_time(ii)) & " sam: " _
            '& Str(Data_Volts(indexit, In_Channel)) & " ref: " & Str(Data_Volts(indexit, Ref_channel)))
            Gain_Volts(0, primary_gain_index, ii + 1) = Data_Volts(indexit, In_Channel(Current_Protocol))
            Gain_Volts(1, primary_gain_index, ii + 1) = Data_Volts(indexit, Ref_Channel_Number(Current_Protocol))

        Next


    End Function



    '    Exit Function

    ''@@@
    '    For pp = 0 To points - 1 Step Number_columns
    '        Column_N = 0

    '        For p = pp To pp + Number_columns - 1
    '            Column_N = Column_N + 1
    ''***************************WITH EXTRA CHANNEL************************

    ''Print(1, Time_Offset + data_time(p) & Chr(9) & Data_Volts(p, In_Channel) & Chr(9))
    '            Gain_Volts_Abs(Column_N) = Data_Volts(p, In_Channel) * (1 - exclude_from_gain And 1)
    '            Gain_Volts_Ref(Column_N) = Data_Volts(p, Ref_channel) * (1 - (exclude_from_gain And 2) / 2)


    '        Next p
    '    Next pp

    '    For ii = 1 To Number_columns

    '        List1.Items.Add("V sample: " & Gain_Volts_Abs(ii) & " ref: " & Gain_Volts_Ref(ii))
    '        If Math.Abs(Gain_Volts_Abs(ii)) < gain_slop And Math.Abs(Gain_Volts_Ref(ii)) < gain_slop Then
    '            List1.Items.Add(" OK ")
    '            Upper_Value = ii
    '        End If
    '    Next



    '    List1.Items.Add(" Upper Value: " & Upper_Value & " TIME = " & Measuring_Pulse_Duration(0, Upper_Value))
    '' figure out the numerical value of the measuring_pulse_duration 




    '    If Upper_Value > 1 Then 'we have ehough values to interpolate the best value
    ''assume that the response is linear between the Upper_Value and the next smallest
    '' time1
    ''calculate slope
    '        dur2 = Duration_Value(Measuring_Pulse_Duration, Upper_Value)
    '        dur1 = Duration_Value(Measuring_Pulse_Duration, Upper_Value - 1)

    '        m = (dur2 - dur1) / (Math.Abs(Gain_Volts_Abs(Upper_Value)) - Math.Abs(Gain_Volts_Abs(Upper_Value - 1)))

    '        List1.Items.Add("slope = " & m & " s/V")
    ''calculate the Y-intercept
    '        b = dur1 / (m - Gain_Volts_Abs(Upper_Value - 1))
    '        List1.Items.Add("interecept = " & b)
    '        target = m * gain_slop '- b
    '        List1.Items.Add("target voltage = " & gain_slop & " --> target duration = " & target & " s")
    '        If target > max_measuring_pulse_duration Then
    '            target = max_measuring_pulse_duration
    '            List1.Items.Add("duration exceeds max allowed, so new target duration = " & target & " s")
    '        End If
    '    Else 'use the highest duration
    '        target = Duration_Value(Measuring_Pulse_Duration, Upper_Value)
    '    End If
    ''make the string to use in program

    '    target_string = Trim(Str$(Math.Round(target * 1000000.0)) & "u")
    '    List1.Items.Add("output string = " & target_string)
    ''Auto_Gain_Calc = target_string

    'End Function
    ' Private Function Duration_Value(ByVal Measuring_Pulse_Duration(,) As String, ByVal Value_Number As Integer) As Single
    '     Dim factor As Single
    '     Dim trimmed As String

    '     trimmed = Trim(Measuring_Pulse_Duration(0, Value_Number))
    '     If VB.Right(trimmed, 1) = "u" Then
    '         factor = 0.000001
    '     ElseIf VB.Right(trimmed, 1) = "m" Then
    '         factor = 0.001
    '     ElseIf VB.Right(trimmed, 1) = "n" Then
    '         factor = 0.000000001
    '     Else
    '         factor = 1
    '     End If
    '     Duration_Value = Val(trimmed) * factor
    ' End Function
    Private Sub subtract_traces(ByRef File_Name_1 As String, ByRef File_Name_2 As String, ByRef File_Name_3 As String)

        Dim Accumulated_Data_X() As Single
        Dim Accumulated_Data_Y() As Single
        Dim Accumulated_Data_Ref() As Single
        Dim data1, data2 As Single
        Dim p, pp, i As Integer
        ReDim Accumulated_Data_X(points)
        ReDim Accumulated_Data_Y(points)
        ReDim Accumulated_Data_Ref(points)
        Dim Reference_Data_Y(points) As Single
        Dim j As String
        Dim temp As Single
        Dim Number_Measuring_Lights As Integer

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)

        'Threading.Thread.sleep(10)

        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i

        Number_Measuring_Lights = Number_Measuring_Lights / 4

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0
        temp = 0
        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        points = i

        'ReDim Baseline_Y(Number_columns)

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data
        FileOpen(2, File_Name_2, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data
        FileOpen(3, File_Name_3, OpenMode.Output, OpenAccess.Default, OpenShare.Shared) 'export subtrated data




        For pp = 0 To points - 1 Step Number_Measuring_Lights

            For p = pp To pp + Number_Measuring_Lights - 1
                Input(2, data2)
                Input(1, data1)
                Print(3, data1 & Chr(9))

                Input(2, data2)
                Input(1, data1)
                Print(3, (data1 - data2) & Chr(9))

                Input(2, data2)
                Input(1, data1)
                Print(3, (data1 - data2) & Chr(9))

                Input(2, data2)
                Input(1, data1)
                Print(3, (data1 - data2) & Chr(9))

            Next p
            PrintLine(3, "")
        Next pp
        FileClose(1)
        'Threading.Thread.sleep(10)

        FileClose(2)
        'Threading.Thread.sleep(10)

        FileClose(3)
        'Threading.Thread.sleep(10)


    End Sub
    Private Sub dice_average(ByRef File_Name_1 As String, ByRef File_Name_2 As String, ByVal dice_length As Integer)
        'list2.Items.Add(File_Name_1 & " -> " & File_Name_2)


        Dim Undiced_Data_X() As Single
        Dim Undiced_Data_Y_raw() As Single
        Dim Undiced_Data_Y_ref() As Single
        Dim Undiced_Data_Y_delta() As Single

        Dim Diced_Data_X() As Single
        Dim Diced_Data_Y_raw() As Single
        Dim Diced_Data_Y_ref() As Single
        Dim Diced_Data_Y_delta() As Single

        'Dim data1, data2 As Single
        Dim p, pp, i, ii, iii As Integer
        Dim j As String
        Dim temp As Single
        Dim Number_Measuring_Lights As Integer
        Dim number_diced As Single

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i
        Number_Measuring_Lights = Number_Measuring_Lights / 4

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0
        temp = 0
        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        points = i

        number_diced = (i / Number_Measuring_Lights) / dice_length

        ReDim Undiced_Data_X(points)
        ReDim Undiced_Data_Y_raw(points)
        ReDim Undiced_Data_Y_ref(points)
        ReDim Undiced_Data_Y_delta(points)

        ReDim Diced_Data_X(points)
        ReDim Diced_Data_Y_raw(points)
        ReDim Diced_Data_Y_ref(points)
        ReDim Diced_Data_Y_delta(points)

        For i = 1 To points
            Diced_Data_X(i) = 0
            Diced_Data_Y_raw(i) = 0
            Diced_Data_Y_ref(i) = 0
            Diced_Data_Y_delta(i) = 0

            Undiced_Data_X(i) = 0
            Undiced_Data_Y_raw(i) = 0
            Undiced_Data_Y_ref(i) = 0
            Undiced_Data_Y_delta(i) = 0

        Next
        'ReDim Baseline_Y(Number_columns)

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data


        i = 0

        For pp = 0 To points - 1 Step Number_Measuring_Lights

            If i = (dice_length) Then i = 0
            iii = 0
            For p = pp To (pp + Number_Measuring_Lights - 1)
                Input(1, Undiced_Data_X(p))
                Input(1, Undiced_Data_Y_raw(p))
                Input(1, Undiced_Data_Y_ref(p))
                Input(1, Undiced_Data_Y_delta(p))
                ii = i * Number_Measuring_Lights + iii
                Diced_Data_X(ii) = Undiced_Data_X(p)
                Diced_Data_Y_raw(ii) = Diced_Data_Y_raw(ii) + (Undiced_Data_Y_raw(p) / number_diced)
                Diced_Data_Y_ref(ii) = Diced_Data_Y_ref(ii) + (Undiced_Data_Y_ref(p) / number_diced)
                Diced_Data_Y_delta(ii) = Diced_Data_Y_delta(ii) + (Undiced_Data_Y_delta(p) / number_diced)
                iii = iii + 1
            Next p
            i = i + 1
        Next pp
        FileClose(1)
        'Threading.Thread.sleep(10)

        FileOpen(3, File_Name_2, OpenMode.Output, OpenAccess.Default, OpenShare.Shared) 'export chopped data


        For pp = 0 To ((dice_length - 1) * Number_Measuring_Lights) Step Number_Measuring_Lights
            For p = pp To (pp + Number_Measuring_Lights - 1)
                Print(3, Diced_Data_X(p) & Chr(9))
                Print(3, Diced_Data_Y_raw(p) & Chr(9))
                Print(3, Diced_Data_Y_ref(p) & Chr(9))
                Print(3, Diced_Data_Y_delta(p))
            Next p
            PrintLine(3, "")

        Next pp
        FileClose(3)
        'Threading.Thread.sleep(10)



    End Sub
    Private Function extract_point(ByRef File_Name_1 As String, ByRef Data_Point As Integer) As Single

        Dim Accumulated_Data_X() As Single
        Dim Accumulated_Data_Y() As Single
        Dim Accumulated_Data_Ref() As Single
        Dim data1, data2 As Single
        Dim p, pp, i As Integer
        ReDim Accumulated_Data_X(points)
        ReDim Accumulated_Data_Y(points)
        ReDim Accumulated_Data_Ref(points)
        Dim Reference_Data_Y(points) As Single
        Dim j As String
        Dim temp As Single
        Dim Number_Measuring_Lights As Integer

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 1
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i
        Number_Measuring_Lights = Number_Measuring_Lights / 4

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0
        temp = 0
        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        points = i

        'ReDim Baseline_Y(Number_columns)

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data

        i = 0

        For pp = 0 To points - 1 Step Number_Measuring_Lights

            For p = pp To pp + Number_Measuring_Lights - 1
                Input(1, data1)
                If i = (((Data_Point + 1) * 4) - 1) Then
                    data2 = data1
                End If
            Next p
            i = i + 1

        Next pp
        FileClose(1)
        'Threading.Thread.sleep(10)

        Return (data2)
    End Function
    Private Function extract_point_raw(ByRef File_Name_1 As String, ByRef Data_Point As Integer) As Single

        Dim Accumulated_Data_X() As Single
        Dim Accumulated_Data_Y() As Single
        Dim Accumulated_Data_Ref() As Single
        Dim data1, data2 As Single
        Dim p, pp, i As Integer
        ReDim Accumulated_Data_X(points)
        ReDim Accumulated_Data_Y(points)
        ReDim Accumulated_Data_Ref(points)
        Dim Reference_Data_Y(points) As Single
        Dim j As String
        Dim temp As Single
        Dim Number_Measuring_Lights As Integer

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 0
        For i = 1 To Len(j)

            'List1.AddItem Mid$(j$, i%, 1) & Asc(Mid$(j$, i%, 1))

            If Asc(Mid(j, i, 1)) = 9 Then
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i
        Number_Measuring_Lights = Number_Measuring_Lights / 4

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)

        i = 0
        temp = 0
        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While
        FileClose(1)
        'Threading.Thread.sleep(10)

        points = i

        'ReDim Baseline_Y(Number_columns)

        FileOpen(1, File_Name_1, OpenMode.Input, OpenAccess.Default, OpenShare.Shared) 'input already saved data

        i = 0

        For pp = 0 To points - 1 Step Number_Measuring_Lights

            For p = pp To pp + Number_Measuring_Lights - 1
                Input(1, data1)
                If i = (((Data_Point + 1) * 4) - 3) Then
                    data2 = data1
                End If
            Next p
            i = i + 1

        Next pp
        FileClose(1)
        'Threading.Thread.sleep(10)

        Return (data2)
    End Function

    Sub Shut_Down()

        cblaster.ChangeCurrentBitmask(0)

    End Sub
    Sub Out_Intensity(ByRef Actinic_Intensity As Short, ByVal blue_actinic As Integer, ByVal far_red As Integer, ByVal saturating_pulse As Integer)

        cblaster.ChangeCurrentBitmask(cblaster.makeBitmask(Actinic_Intensity, 0, 0, 0, 0, far_red, 0, blue_actinic, 0, 0))


        ProgressBar3.Value = Actinic_Intensity
        actinic_label.Text = Actinic_Intensity
    End Sub

    Public Sub halt_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles halt.Click

        Halt_Script = True
        Call Shut_Down()

    End Sub

    Public Sub LEDS_Off_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles LEDS_Off.Click
        Call Shut_Down()
        ProgressBar3.Value = 0
        actinic_label.Text = 0
    End Sub

    Public Sub Load_Script_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles Load_Script.Click
        Call Do_Load_Script()
    End Sub
    Private Sub Do_Load_Script()
        CommonDialog1Open.FileName = "*.txt"
        'CommonDialog1Save.FileName = "*.txt"
        'CommonDialog1Open.Filter = "Script Files|*.scr|Text Files|*.txt"
        CommonDialog1Open.Filter = "Text Files|*.txt"
        CommonDialog1Save.Title = "OPEN SCRIPT"

        CommonDialog1Save.FileName = CommonDialog1Open.FileName
        CommonDialog1Open.ShowDialog()

        'CommonDialog1Save.Filter = "Script Files|*.scr|Text Files|*.txt"
        script_file = CommonDialog1Open.FileName
        Script_Label.Text = script_file

    End Sub

    Public Sub Run_DMK_Script_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles Run_DMK_Script.Click
        run_another_script = True
        'script_file = Script_Label.Text
        'Me.Text = "Kinetic Spectrophotometer: RUNNING"
        Jabber_jabber.FilesListBox.Items.Add("running script: " & script_file)
        Me.BackColor = Color.Green
        list_of_base_additions.Clear()

        While run_another_script = True
            Call run_the_script()
        End While
        'Me.Text = "Kinetic Spectrophotometer: IDLE"
        Me.BackColor = Color.Gray
    End Sub

    ' pass false to this function to make it run int test mode
    Private Sub run_the_script()
        run_the_script(True)
    End Sub
    Private Sub run_the_script(ByVal ForReal As Boolean)
        ' if ForReal is true, then the script will run normally. If ForReal is false, then the script will run in test mode
        'On Error GoTo err_it
        run_another_script = False
        If Trace_Running = True Then Exit Sub

        'this sub parses the script into arrays and feeds the arrays into the Multi_Trace sub

        'Dim Max_Loops As Integer
        'Dim script_file As String
        '        
        'Dim Script() As String
        Dim temp As String = ""
        Dim i, ii As Short


        Const Max_Loops As Short = 100


        Dim L_Beg(Max_Loops) As Object
        Dim L_End(Max_Loops) As Object
        Dim L_Number(Max_Loops) As Object
        Dim L_counter(Max_Loops) As Short
        Dim Total_Script_Count As Short
        Dim Loop_Number As Short
        Dim Wait_Time As Single

        ' dimentions involving DIRK traces


        Dim Measuring_Pulse_Duration(,) As String = {}
        'Dim Number_Dirk_Repeats() As Short = {0, 0, 0}
        'Dim Dirk_Measuring_Interval() As String = {"0", "0", "0"} ' time between measuring pulses
        Dim Baseline_Measuring_Interval() As String = {"0", "0", "0"}
        Dim Recovery_Measuring_Interval() As String = {"0", "0", "0"}
        Dim Baseline_Number_Points() As Short = {0, 0, 0}
        Dim Recovery_Number_Points() As Short = {0, 0, 0}
        Dim Dirk_Number_Points() As Short = {0, 0, 0}
        'Dim Pre_Delay_Time As String

        Dim trace1, trace2, trace3 As Integer
        Dim simmer_number As Integer
        Dim simmer_q_switch As Integer


        Dim q_switch_time As Integer

        Dim In_Channel() As Short = {0, 0, 0}
        Dim Aux_In_Channel() As Short = {0, 0, 0}
        Dim Ref_Channel_Number() As Short = {0, 0, 0}
        Dim Aux_Ref_Channel() As Short = {0, 0, 0}

        'Dim Baseline_Actinic_Intensity() As Short = {0, 0, 0}
        'Dim Recovery_Actinic_Intensity() As Short = {0, 0, 0}
        'Dim Dirk_Actinic_Intensity() As Short = {0, 0, 0}

        'Dim Dirk_Data_X() As Single
        'Dim Dirk_Data_Y() As Single
        'Dim Dirk_Points As Short
        Dim Time_Mode() As String = {"0", "0", "0"}
        Dim trace_note() As String = {"0", "0", "0"}
        Dim trace_label() As String = {"0", "0", "0"}

        Dim ulstat As Short
        Dim temp_string As String

        Dim Number_Protocols As Short
        Dim Save_Mode() As Short = {0, 0, 0}
        Dim file_name_sequence() As short = {0,0,0}

        Dim Times_Averaged() As Short = {0, 0, 0}
        'Dim Fluorescence_Times_Averaged() As Short
        Dim Wavelength() As String = {"none", "none", "none"}
        Dim Gain() As Integer = {0, 0, 0}
        Dim Blocking_Filter() As String = {"0", "0", "0"}
        Dim Give_Flash() As Short = {0, 0, 0}
        Dim Wheel_Position() As Short = {0, 0, 0}
        Dim Detector_Gain() As Short = {0, 0, 0}

        Dim File_Name() As String = {"temp", "temp", "temp"}
        Dim Aux_File_Name() As String = {"temp", "temp", "temp"}

        Dim Current_Protocol As Short

        Dim wt As Short

        Dim Number_Traces As Short

        'Dim set_gain_for_this_measuring_light As String
        Dim Current_Intensity As Short = 0
        Dim Far_Red As Short
        Dim Blue_Actinic As Short

        ' END dimentions involving DIRK traces

        ' start DIm with MULTI-traces
        'Call Multi_Trae(Number_Loops, Intensity(), Number_Pulses(), Number_Measuring_Lights, Measuring_Light(), Measuring_Interval(), Primary_Gain(), Measuring_Pulse_Duration)
        Dim alternate_reference_file_name() As String = {"", "", ""}
        Dim M_Number_Loops() As Short = {0, 0, 0}
        Dim M_Intensity(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim M_Far_Red(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim M_Blue_Actinic(,) As Short = {{0, 0, 0}, {0, 0, 0}}

        Dim protocol_label() As String = {"0", "0", "0"}

        Dim Pre_Pulse(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim Pre_Delay(,) As String = {{"0", "0", "0"}, {"0", "0", "0"}}
        Dim Pre_Pulse_Time(,) As String = {{"0", "0", "0"}, {"0", "0", "0"}}
        Dim M_Number_Pulses(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim M_Number_Measuring_Lights() As Short = {0, 0, 0}
        Dim M_Measuring_Light(,) As Short = {{0, 0, 0}, {0, 0, 0}}

        Dim M_Measuring_Light_Names(,) As String = {{"1", "2", "3"}, {"4", "5", "6"}}


        'Dim M_Wasp(,) As Short = {{0, 0, 0}, {0, 0, 0}}

        Dim m_plot_specific_traces(,) As Short ' = {{0, 0}, {0, 0}}
        Dim M_Take_Data(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim M_choke_actinic(,) As Short = {{0, 0, 0}, {0, 0, 0}}

        Dim M_Measuring_Interval(,) As String = {{0, 0, 0}, {0, 0, 0}}
        Dim L_Measuring_Interval(,,) As String = {{{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}}
        'Dim M_Primary_Gain(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim Xe_Flash(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim q_switch(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        'the following line is deleted in the Exxon version
        'Dim M_Reference_Gain(,) As Short = {{0, 0, 0}, {0, 0, 0}}
        Dim Baseline_Start() As Short = {0, 0, 0}
        Dim Baseline_End() As Short = {0, 0, 0}
        Dim iii As Short
        Dim M_Temp As Short
        Dim Max_Number_Loops As Short

        Max_Number_Loops = 100

        If ForReal Then
            Jabber_jabber.List1.Items.Clear()
        Else
            Jabber_jabber.List1.Items.Add("reading script file in test mode...")
            Jabber_jabber.List1.TopIndex = List1.Items.Count - 1
        End If

        'Measuring_Pulse_Duration
        'Intensity()
        'Dim Max_Volts As Single
        'Dim gain_temp As Single

        'Dim gain_set_temp As Integer
        Dim Start_Time As Single
        Dim Zero_Time As Single
        Dim gain_data(10) As Single
        Dim Array_Value(10, 100) As Object
        Dim Array_Number As Short
        Dim Array_Max_Index(10) As Object
        Dim Array_Index(10) As Short
        Dim points, temp_s As Integer
        Dim Last_Servo As String = ""
        Dim iv As Short
        Dim smooth() As Short
        'Dim offset As Single

        Dim Wheel(25) As Single 'contains the positions of the wheel openings for the servo
        Dim word_hold As String

        Dim exclude_from_gain As Integer
        'FileOpen(3, "c:\wheel_position.dat", OpenMode.Input)
        'For iii = 0 To 25
        ' Input(3, Wheel(iii))
        ' Next iii
        'FileClose(3)

        Zero_Time = VB.Timer() 'indicates the start of the trace
        Halt_Script = False
        'Dim auto_gain_value As String


        'File_Name = "c:\test.dat"
        'Measuring_Pulse_Duration = "12u"
        'Baseline_Measuring_Interval = "100u"
        'Dirk_Measuring_Interval = "100u"
        'Recovery_Measuring_Interval = "1000u"
        'Dirk_Number_Points = 100
        'Baseline_Number_Points = 100
        'Recovery_Number_Points = 500
        'Actinic_Intensity = 0
        'Dirk_Actinic_Intensity = 16

        'Dim Number_Fluorescence_Traces As Integer


        'Dim Induction_Measuring_Pulse_Duration()
        'Dim Induction_Baseline_Measuring_Interval()
        'Dim Induction_Baseline_Measuring_Points()
        'Dim Induction_Measuring_Interval()

        'Dim Induction_Number_Points() As Integer


        'Dim Initial_Actinic_Intensity() As Integer

        'Dim Final_Acinic_Intensity() As Integer

        'Dim Induction_Gain() As Integer

        'Dim Measuring_Mode() As String = {"520", "820", "fluor"}

        'fluor
        '520
        '820
        'If script_file = "" Then
        '    Call Do_Load_Script
        'End If
        ' see how long the script is

        'File_Name = "c:\test.dat"

        'script_file = Script_Label.Text


        If script_file = "" Then
            CommonDialog1Open.Title = "LOAD SCRIPT FILE"
            CommonDialog1Save.Title = "LOAD SCRIPT FILE"

            CommonDialog1Open.FileName = "*.txt"
            'CommonDialog1Save.FileName = "*.txt"
            'Else
            'CommonDialog1Open.FileName = script_file
            'End If

            CommonDialog1Open.ShowDialog()
            'CommonDialog1Save.FileName = CommonDialog1Open.FileName
            'UPGRADE_WARNING: Filter has a new behavior. Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="9B7D5ADD-D8FE-4819-A36C-6DEDAF088CC7"'
            'CommonDialog1Open.Filter = "Script Files|*.scr|Text Files|*.txt"
            'CommonDialog1Save.Filter = "Script Files|*.scr|Text Files|*.txt"
            CommonDialog1Open.Filter = "Text Files|*.txt"
            CommonDialog1Save.Filter = "Text Files|*.txt"

            script_file = CommonDialog1Open.FileName
            Script_Label.Text = script_file

            If CommonDialog1Open.CheckFileExists = False Then
                System.Windows.Forms.MessageBox.Show("no such file")
                script_file = ""
                Exit Sub
            End If
            If script_file = "*.*" Or script_file = "*.txt" Then
                script_file = ""
                Exit Sub
            End If

        End If

        FileOpen(1, script_file, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        FileOpen(3, rootdata_directory & "include_temp.txt", OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
        Total_Script_Count = 0
        While EOF(1) = False
            Input(1, temp)
            temp = LCase(temp)
            If VB.Left(temp, 1) = ">" Then  'merge another file into this one
                temp = VB.Mid(temp, 2, (VB.Len(temp) - 1))
                FileOpen(2, temp, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
                While EOF(2) = False
                    Input(2, temp)
                    If LCase(temp) = ">" Then  'merge another file into this one
                        System.Windows.Forms.MessageBox.Show("Error: The program cannot handle recursive INCLUDE (>) statements", "ERROR", MessageBoxButtons.OK, MessageBoxIcon.Error)
                        FileClose(1)
                        FileClose(2)
                        FileClose(3)
                        Exit Sub
                    End If
                    PrintLine(3, temp)
                End While
                FileClose(2)
            Else
                PrintLine(3, temp)
            End If
        End While
        FileClose(1)
        FileClose(3)

        FileOpen(1, rootdata_directory & "include_temp.txt", OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        FileOpen(3, rootdata_directory & "parse_temp.txt", OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

        Total_Script_Count = 0
        While EOF(1) = False
            Input(1, temp)
            temp = Trim(temp)
            temp = LCase(temp)
            If temp = "" Then GoTo skipit
            If VB.Left(temp, 1) = "'" Then GoTo skipit
            word_hold = ""
            i = 1
            While i < VB.Len(temp) + 1  'in this section, we parse the script into nuggets, separated by CRLF, like the old script files

                If VB.Mid(temp, i, 1) = "," Or VB.Mid(temp, i, 1) = "(" Or VB.Mid(temp, i, 1) = ")" Or VB.Mid(temp, i, 1) = ";" Then  'ignore parentheses or commas
                    PrintLine(3, word_hold)
                    word_hold = ""
                    'there is a problem here! entering a negative number causes the number to be split!
                    ElseIf VB.Mid(temp, i, 1) = "+" Or VB.Mid(temp, i, 2) = "--" Or VB.Mid(temp, i, 1) = "*" Or VB.Mid(temp, i, 1) = "=" Or VB.Mid(temp, i, 1) = "/" Then  'this means we are doing math
                       PrintLine(3, word_hold)
                       PrintLine(3, VB.Mid(temp, i, 1))
                       i = i + 1
                       word_hold = VB.Mid(temp, i, 1)
                ElseIf VB.Mid(temp, i, 1) = "'" Then  'this means the rest of the line is a comment
                    GoTo skip_comment
                Else
                    word_hold = word_hold + VB.Mid(temp, i, 1)

                End If

                i = i + 1
            End While
skip_comment:
            If word_hold <> "" Then PrintLine(3, word_hold)
            'If VB.Left(temp, 1) = "@" Then PrintLine(3, "0") 'add another line to hold the value of a variable in the script
skipit:
        End While

        FileClose(3)
        FileClose(1)

        'Exit Sub

        '        temp = Trim(temp)
        '        If LCase(temp) = "include" Then  'merge another file into this one
        'Input(1, temp)
        'FileOpen(2, temp, OpenMode.Input)
        'While EOF(2) = False
        ' Input(2, temp)
        ' temp = Trim(temp)
        ' If VB.Left(temp, 1) <> "'" And temp <> "" Then ' ' designates a remark
        ' Total_Script_Count = Total_Script_Count + 1
        ' End If
        ' End While
        ' FileClose(2)
        ' End If
        ' If VB.Left(temp, 1) <> "'" And temp <> "" Then ' ' designates a remark
        ' Total_Script_Count = Total_Script_Count + 1
        ' End If
        '        End While
        ' FileClose(1)
        ' FileClose(3)

        'Exit Sub

        FileOpen(1, rootdata_directory & "parse_temp.txt", OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        ' FileOpen(3, "c:\temp.txt", OpenMode.Output)

        Total_Script_Count = 0
        While EOF(1) = False
            Input(1, temp)
            temp = Trim(temp)
            If VB.Left(temp, 1) <> "'" And temp <> "" Then ' ' designates a remark
                Total_Script_Count = Total_Script_Count + 1
            End If
        End While
        FileClose(1)


        ReDim Script(Total_Script_Count)

        ' now input script into the array "SCRIPT"

        FileOpen(1, rootdata_directory & "parse_temp.txt", OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        Script_Counter = 0
        While EOF(1) = False
            Input(1, temp)
            temp = Trim(temp)
            If VB.Left(temp, 1) <> "'" And temp <> "" Then ' ' designates a remark or empty line
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1
                Script(Script_Counter) = LCase(temp)
                If LCase(temp) = "include" Then  'include another file into this one
                    Input(1, temp)
                    FileOpen(2, temp, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
                    While EOF(2) = False
                        Input(2, temp)
                        temp = Trim(temp)
                        If VB.Left(temp, 1) <> "'" And temp <> "" Then ' ' designates a remark or empty line
                            Call Advance_Script()
                            'Script_Counter = Script_Counter + 1
                            Script(Script_Counter) = LCase(temp)
                        End If
                    End While
                    FileClose(2)
                End If
            End If
        End While
        FileClose(1)

        'Now run the script

        Loop_Number = 0
        Script_Counter = 0
        Sub_Return_Level = 0

        While Script_Counter < Total_Script_Count And Halt_Script = False
            Call Advance_Script()
            'Script_Counter = Script_Counter + 1
            'List1.AddItem(">" & Script(Script_Counter))

            If C_Script = "lb" Then 'designates the beginning of a loop
                ' this script command has the following form:
                'lb
                '1 (loop number, max of 20)
                '10 (number of loops to do)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Loop_Number = Val(C_Script) 'set number of loops to do
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                L_Number(Loop_Number) = Val(C_Script) 'loop back at the next script_counter

                L_Beg(Loop_Number) = Script_Counter '+ 1 'loop back at the next script_counter

                '            End If

            ElseIf C_Script = "le" Then 'designates the beginning of a loop
                ' this script command has the following form:
                'lb
                '1 (loop number, max of 20)
                '10 (number of loops to do)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Loop_Number = Val(C_Script) 'set number of loops to do

                L_Number(Loop_Number) = L_Number(Loop_Number) - 1 'decrement L_number()

                If ForReal Then Jabber_jabber.List1.Items.Add("LOOP:" & Loop_Number & "TO GO:" & L_Number(Loop_Number))
                'Scrolling text box fix: make the box jump to the bottom after adding an item
                List1.TopIndex = List1.Items.Count - 1
                If L_Number(Loop_Number) > 0 Then
                    Script_Counter = L_Beg(Loop_Number) 'loop back at the next script_counter
                End If
                'End If



            ElseIf C_Script = "sub" Then 'designates a call to a subroutine
                ' this script command has the following form:
                'sub, name_of_subroutine

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                sub_name = C_Script & "|" 'set number of loops to do
                Sub_Return_Level = Sub_Return_Level + 1
                Sub_Return_Pointer(Sub_Return_Level) = Script_Counter + 1


                'Label7.Text = ("not_found:" & sub_name)

                ' find the subroutine start point
                For iv = 1 To Total_Script_Count

                    If Script(iv) = sub_name Then
                        Script_Counter = iv
                        '       Label7.Text = ("found: " & sub_name)

                    End If
                Next iv

                'End If

            ElseIf C_Script = "return" Then 'designates the beginning of a loop
                Script_Counter = Sub_Return_Pointer(Sub_Return_Level) - 1
                Sub_Return_Level = Sub_Return_Level - 1
                'End If

            ElseIf C_Script = "end" Then 'designates the end of the script
                'Label7.Text = "end"
                Exit Sub
                'End If

            ElseIf VB.Left(C_Script, 1) = "!" Then 'designates that we are defining a variable

                'Script_Counter = Script_Counter + 1 'increment the script counter

                '                Script_Counter = Script_Counter + 1 'increment the script counter
                'If VB.Mid(Script(Script_Counter), 2, 1) <> "@" Then
                ' System.Windows.Forms.MessageBox.Show("Variables must begin with '@'", "ERROR", MessageBoxButtons.OK, MessageBoxIcon.Error)
                ' Exit Sub
                'End If
                number_variables = number_variables + 1
                'iv = VB.Len(Script(Script_Counter)) - 1
                user_variable_name(number_variables) = VB.Mid(C_Script, 2, (VB.Len(C_Script) - 1))
                user_variable_index(number_variables) = 1
                If ForReal Then List1.Items.Add(user_variable_name(number_variables))
                'Scrolling text box fix: make the box jump to the bottom after adding an item
                List1.TopIndex = List1.Items.Count - 1
                'Call Advance_Script()
                'Script_Counter = Script_Counter + 1
                'i = Val(Script(Script_Counter))
                'List1.Items.Add(user_variable_name(number_variables))
                'End If

            ElseIf VB.Left(C_Script, 1) = "#" Then 'designates that we are doing arithmatic on a variable
                Dim j As Single
                Dim k As Integer
                Dim var As String
                var = VB.Right(C_Script, (Len(C_Script) - 1))
                'Label7.Text = var
                j = 0
                For i = 1 To number_variables
                    If var = user_variable_name(i) Then
                        j = i
                    End If
                Next i
                If j = 0 Then
                    System.Windows.Forms.MessageBox.Show("Error: I cannot find that variable", "ERROR", MessageBoxButtons.OK, MessageBoxIcon.Error)
                    Exit Sub
                End If

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                k = Val(C_Script)
                user_variable_index(j) = k

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                Select Case C_Script
                    Case "+"
                        If ForReal Then List1.Items.Add("+")
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        user_variable(j, k) = Val(user_variable(j, k)) + Val(C_Script)
                    Case "-"
                        If ForReal Then List1.Items.Add("-")
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        user_variable(j, k) = user_variable(j, k) - Val(C_Script)

                    Case "*"
                        If ForReal Then List1.Items.Add("*")
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        user_variable(j, k) = Val(user_variable(j, k)) * Val(C_Script)

                    Case "/"
                        If ForReal Then List1.Items.Add("/")
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        user_variable(j, k) = Val(user_variable(j, k)) / Val(C_Script)

                    Case "="
                        If ForReal Then List1.Items.Add("=")
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        user_variable(j, k) = Val(C_Script)

                    Case Else
                        System.Windows.Forms.MessageBox.Show("Error: invalid operator", "ERROR", MessageBoxButtons.OK, MessageBoxIcon.Error)
                        Exit Sub
                End Select
                If ForReal Then List1.Items.Add(user_variable_name(j) & ":" & k & " = " & user_variable(j, k))
                'Scrolling text box fix: make the box jump to the bottom after adding an item
                List1.TopIndex = List1.Items.Count - 1
                'End If
            ElseIf C_Script = "wait" Then 'designates the beginning of a loop
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If C_Script = "v" Then
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Array_Number = Val(C_Script)

                    Wait_Time = Val(Array_Value(Array_Number, Array_Index(Array_Number)))
                Else
                    Wait_Time = Val(C_Script) 'set wait time
                End If
                If ForReal Then Call Hold_on_There(Wait_Time)
            ElseIf C_Script = "gain_slop" Or C_Script = "gain_target" Then
                Call Advance_Script()
                gain_slop = Val(C_Script) '

            ElseIf C_Script = "ref_channel" Or C_Script = "ref_mode" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                If Val(C_Script) > 0 Then
                    ref_mode = 1
                    Jabber_jabber.List1.Items.Add("ref_mode set to 1: use ref channel data to calculate deltaA")

                Else
                    ref_mode = 0
                    Jabber_jabber.List1.Items.Add("ref_mode set to 0: use baseline data to calculate deltaA")
                End If


            ElseIf C_Script = "set_base_file" Then
                If ForReal Then Call set_base_file_name()


            ElseIf C_Script = "auto_base_file" Then
                If ForReal Then

                    Base_File_Name = Trim("IDEA" & DateString & TimeOfDay)
                    For i = 1 To Base_File_Name.Length
                        If Mid(Base_File_Name, i, 1) = ":" Then
                            Mid(Base_File_Name, i, 1) = "_"
                        End If
                    Next

                    Base_File_Name = Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) & "\Documents\idea_spec_test_data\" & Base_File_Name



                    If Directory.Exists(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) & "\Documents\idea_spec_test_data") = False Then
                        Directory.CreateDirectory(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) & "\Documents\idea_spec_test_data")
                    End If

                    List1.Items.Add(Base_File_Name)

                End If


            ElseIf C_Script = "set_base_file_pre" Then
                If ForReal Then

                    Call Advance_Script()
                    Base_File_Name = C_Script
                End If


            ElseIf C_Script = "note_query" Then
                If ForReal Then

                    Note_File_Name = Trim(Base_File_Name & "_notes.txt")
                    take_a_note.Show()
                    wait_for_notes = True
                    While wait_for_notes = True
                        System.Windows.Forms.Application.DoEvents()

                    End While
                End If
            ElseIf C_Script = "record_events" Then  'sets whether to adds events to the note file
                Call Advance_Script()
                If ForReal Then
                    If C_Script > 0 Then
                        take_notes = True
                    Else
                        take_notes = False
                    End If
                End If
            ElseIf C_Script = "record_script" Then  'adds text to the note file
                If ForReal Then
                    note_text = note_text & "*******************************************************************" & ControlChars.CrLf
                    note_text = note_text & "TEXT OF SCRIPT FILE USED:" & ControlChars.CrLf

                    FileOpen(1, script_file, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
                    While EOF(1) = False
                        temp = LineInput(1)
                        note_text = note_text & temp & ControlChars.CrLf

                    End While

                    note_text = note_text & ControlChars.CrLf & "*******************************************************************" & ControlChars.CrLf
                    FileClose(1)
                    'Threading.Thread.Sleep(10)

                    FileOpen(7, Note_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
                    Write(7, note_text)
                    FileClose(7)
                    'Threading.Thread.Sleep(10)

                End If
            ElseIf C_Script = "make_note" Then  'adds text to the note file
                If ForReal Then note_text = note_text & ControlChars.CrLf & "*******************************************************************" & ControlChars.CrLf
                Call Advance_Script()
                If ForReal Then
                    note_text = note_text & Date.Now & C_Script & ControlChars.CrLf
                    note_text = note_text & ControlChars.CrLf & "*******************************************************************" & ControlChars.CrLf
                    FileOpen(7, Note_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)

                    Write(7, note_text)
                    FileClose(7)
                    'Threading.Thread.Sleep(10)

                End If

            ElseIf C_Script = "shut_down" Then
                If ForReal Then Call Shut_Down()

            ElseIf C_Script = "append_base" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                i = Current_Trace
                Base_File_Name = Base_Base_File_Name + C_Script
                'List1.Items.Add("**" & Base_File_Name)
                If list_of_base_additions.Contains(C_Script) Then
                    List1.Items.Add("found " & C_Script)
                    File_Name(i) = Trim(Base_File_Name & VB6.Format(i, "0000") & ".dat")
                    Aux_File_Name(i) = Trim(Base_File_Name & "_aux_" & VB6.Format(i, "0000") & ".dat")
                Else

                    For i = 0 To Number_Traces ' put in some default values
                        File_Name(i) = Trim(Base_File_Name & VB6.Format(i, "0000") & ".dat")
                        Aux_File_Name(i) = Trim(Base_File_Name & "_aux_" & VB6.Format(i, "0000") & ".dat")
                        Times_Averaged(i) = 0
                    Next i

                End If
                list_of_base_additions.Add(C_Script) ' add to the list of values used
                For Each bas As String In list_of_base_additions
                    List1.Items.Add(bas)
                Next

            ElseIf C_Script = "link" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                script_file = C_Script

                Trace_Running = False
                run_another_script = True
                Exit Sub


            ElseIf C_Script = "set_lambda" Then
                Call Advance_Script()
                TextBox4.Text = C_Script

                Form3.moveTo.Text = Val(C_Script) 'TextBox4.Text
                Label14.Text = "moving"

                Form3.moveIt()

                While Form3.moving = True
                    System.Windows.Forms.Application.DoEvents()
                End While

                Label14.Text = "done"

                'Script_Counter = Script_Counter + 1 'increment the script counter

                'lambdasettext.Text = C_Script
                'Set_Lambda(Val(C_Script))
            ElseIf C_Script = "set_linear" Then

                Call Advance_Script()
                TextBox4.Text = C_Script

                Form3.linPosTo.Text = Val(C_Script) 'TextBox4.Text
                Label14.Text = "moving"

                Form3.lMoveIt()

                While Form3.moving = True
                    System.Windows.Forms.Application.DoEvents()
                End While

                Label14.Text = "done"

            ElseIf C_Script = "set_linear_index" Then

                Call Advance_Script()
                TextBox4.Text = C_Script

                Form3.linMoveToIndexIndex.Text = Val(C_Script) 'TextBox4.Text
                Label14.Text = "moving"

                Form3.linMoveToIndexf()

                While Form3.moving = True
                    System.Windows.Forms.Application.DoEvents()
                End While

                Label14.Text = "done"

            ElseIf C_Script = "measuring_pulse_duration" Then
                'Dim xxx As String
                For n_light As Integer = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Measuring_Pulse_Duration(Current_Protocol, n_light) = C_Script
                    ' ERROR CATCH: If the pulse duration is greator than 400 us, then the LED will DIE!!!
                    ' The following block will catch this problem
                    ' ///// begin pulse duration catch /////
                    Dim trimmed As String
                    Dim factor As Double
                    Dim duration As Double
                    trimmed = Trim(Measuring_Pulse_Duration(Current_Protocol, n_light))
                    If VB.Right(trimmed, 1) = "u" Then
                        factor = 0.000001
                    ElseIf VB.Right(trimmed, 1) = "m" Then
                        factor = 0.001
                    ElseIf VB.Right(trimmed, 1) = "n" Then
                        factor = 0.000000001
                    Else
                        factor = 1
                    End If
                    duration = Val(trimmed) * factor
                    If duration > (400 * 0.000001) Then
                        Dim keepgoing As Boolean
                        keepgoing = ShowScriptErrorMessage(Script_Counter, Script, "DANGER: measuring pulse length" & Str(duration) & " is too long and will probably cause permanent damage to the measuring LED.")
                        'keepgoing = ShowScriptErrorMessage(Script_Counter, Script, "DANGER: measuring pulse length is too long and will probably cause permanent damage to the measuring LED.")
                        If keepgoing = False Then Exit Sub
                    End If
                    'List1.Items.Add("protocol: " + Str(Current_Protocol) + "%%%:" + Measuring_Pulse_Duration(Current_Protocol, n_light))
                Next n_light
                ' ///// end pulse duration catch /////
                'List1.AddItem "mpi" & Measuring_Pulse_Duration(Current_Protocol)

            ElseIf C_Script = "protocol_label" Then
                Call Advance_Script()
                protocol_label(Current_Protocol) = C_Script


            ElseIf C_Script = "m_number_loops" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                M_Number_Loops(Current_Protocol) = Val(C_Script) '
                For iii = 1 To M_Number_Loops(Current_Protocol)  'set to take data on all loops; can be overridden with m_take_data
                    M_Take_Data(Current_Protocol, iii) = 1
                Next
            ElseIf C_Script = "m_separate_reference_trace" Then
                Call Advance_Script()
                If Val(C_Script) > 0 Then
                    alternate_reference_file_name(Current_Trace) = File_Name(Val(C_Script))
                Else
                    alternate_reference_file_name(Current_Trace) = ""
                End If


            ElseIf C_Script = "baseline_start" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Baseline_Start(Current_Protocol) = Val(C_Script) '

            ElseIf C_Script = "beep" Then
                Dim beep_freq, beep_dur As Integer
                Call Advance_Script()
                beep_freq = Val(C_Script) '
                Call Advance_Script()
                beep_dur = Val(C_Script) '

                'Beep()
                System.Console.Beep(beep_freq, beep_dur)

            ElseIf C_Script = "baseline_end" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Baseline_End(Current_Protocol) = Val(C_Script) '

            ElseIf C_Script = "laser" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                laser = Val(C_Script) '

                'Dim M_Number_Loops As Integer

            ElseIf C_Script = "m_intensity" Then


                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter

                    If C_Script = "v" Then
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        Array_Number = Val(C_Script)
                        M_Intensity(Current_Protocol, M_Temp) = Val(Array_Value(Array_Number, Array_Index(Array_Number)))
                    Else
                        M_Intensity(Current_Protocol, M_Temp) = Val(C_Script) 'set wait time
                    End If
                Next M_Temp

            ElseIf C_Script = "m_far_red" Then

                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Far_Red(Current_Protocol, M_Temp) = Val(C_Script) 'set wait time
                Next M_Temp

            ElseIf C_Script = "m_blue_actinic" Or C_Script = "m_flash" Then


                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Blue_Actinic(Current_Protocol, M_Temp) = Val(C_Script) 'set wait time

                Next M_Temp


            ElseIf C_Script = "m_leaf_actinic" Then


                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Blue_Actinic(Current_Protocol, M_Temp) = Val(C_Script) 'set wait time

                Next M_Temp

                'Dim M_Intensity() As Integer
            ElseIf C_Script = "xe_flash" Then

                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Xe_Flash(Current_Protocol, M_Temp) = CShort(C_Script)
                    If ForReal Then List1.Items.Add("xe_added" & Xe_Flash(Current_Protocol, M_Temp))
                    'Scrolling text box fix: make the box jump to the bottom after adding an item
                    List1.TopIndex = List1.Items.Count - 1
                Next M_Temp

            ElseIf C_Script = "q_switch" Then

                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    q_switch(Current_Protocol, M_Temp) = CShort(C_Script)
                    'List1.Items.Add("xe_added" & Xe_Flash(Current_Protocol, M_Temp))
                Next M_Temp


            ElseIf C_Script = "m_pulse_set" Then

                For M_Temp = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Number_Pulses(Current_Protocol, M_Temp) = Val(C_Script)
                Next M_Temp

                'Dim M_Number_Pulses() As Integer

            ElseIf C_Script = "m_number_measuring_lights" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                M_Number_Measuring_Lights(Current_Protocol) = Val(C_Script)

                'Dim M_Number_Measuring_Lights As Integer
            ElseIf C_Script = "m_take_data" Then
                Dim mm As Integer
                For mm = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()

                    M_Take_Data(Current_Protocol, mm) = Val(C_Script)

                Next mm

            ElseIf C_Script = "m_plot_trace" Then

                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    m_plot_specific_traces(Current_Protocol, M_Temp) = Val(C_Script)
                Next M_Temp

            ElseIf C_Script = "m_measuring_light" Then

                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Measuring_Light(Current_Protocol, M_Temp) = Val(C_Script)
                Next M_Temp

            ElseIf C_Script = "m_measuring_light_names" Then

                For M_Temp = 0 To M_Number_Measuring_Lights(Current_Protocol) - 1
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Measuring_Light_Names(Current_Protocol, M_Temp) = C_Script
                Next M_Temp

            ElseIf C_Script = "wasp" Then


                'For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '    Call Advance_Script()
                '    'Script_Counter = Script_Counter + 1 'increment the script counter
                '    M_Wasp(Current_Protocol, M_Temp) = Val(C_Script) 'set wait time

                'Next M_Temp

            ElseIf C_Script = "m_choke_actinic" Then 'one can now disable the actinic light specifically when running each measuring light; the dark interval will last for the duration of the delay time and measuring pulse etc  

                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_choke_actinic(Current_Protocol, M_Temp) = Val(C_Script)
                Next M_Temp


                'Dim M_Measuring_Light() As Integer

                '            ElseIf c_script = "m_measuring_interval" Then

                '                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                ' Script_Counter = Script_Counter + 1 'increment the script counter
                ' M_Measuring_Interval(Current_Protocol, M_Temp) = (c_script)
                ' Next M_Temp

            ElseIf C_Script = "m_measuring_interval" Or C_Script = "l_measuring_interval" Then
                Dim mm As Integer
                For mm = 1 To M_Number_Loops(Current_Protocol)
                    For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1 'increment the script counter
                        L_Measuring_Interval(Current_Protocol, M_Temp, mm) = (C_Script)
                    Next M_Temp
                Next mm
            ElseIf C_Script = "m_pre_pulse_light" Or C_Script = "pre_pulse_light" Then
                Call Advance_Script()
                pre_pulse_light = Val(C_Script)
            ElseIf C_Script = "m_pre_pulse_measuring_light" Then
                Call Advance_Script()
                If Val(C_Script) > 0 Then
                    pre_pulse_light = pre_pulse_light Or Val(C_Script) * &H10
                Else
                    pre_pulse_light = pre_pulse_light And &HFF0F
                End If
                Jabber_jabber.List1.Items.Add("pre_pulse_light = " + Str(pre_pulse_light))
            ElseIf C_Script = "m_pre_pulse_actinic_intensity" Then
                Call Advance_Script()
                pre_pulse_light = pre_pulse_light And &HFF 'reset the upper byte, then in the next step refill

                pre_pulse_light = pre_pulse_light Or Val(C_Script) * &H100


            ElseIf C_Script = "m_pre_pulse_xe" Then
                Call Advance_Script()
                If Val(C_Script) = 1 Then
                    pre_pulse_light = pre_pulse_light Or 8
                Else
                    pre_pulse_light = pre_pulse_light And (&HFFFF - 8)
                End If
                Jabber_jabber.List1.Items.Add("pre_pulse_light = " + Str(pre_pulse_light))

            ElseIf C_Script = "m_pre_pulse_blue" Then
                Call Advance_Script()
                If Val(C_Script) = 1 Then
                    pre_pulse_light = pre_pulse_light Or 2
                Else
                    pre_pulse_light = pre_pulse_light And (&HFFFF - 2)
                End If
                Jabber_jabber.List1.Items.Add("pre_pulse_light = " + Str(pre_pulse_light))

            ElseIf C_Script = "m_pre_pulse_far_red" Then
                Call Advance_Script()
                If Val(C_Script) = 1 Then
                    pre_pulse_light = pre_pulse_light Or 4
                Else
                    pre_pulse_light = pre_pulse_light And (&HFFFF - 4)
                End If
                Jabber_jabber.List1.Items.Add("pre_pulse_light = " + Str(pre_pulse_light))

                'pre_pulse_light = pre_pulse_light + 4 * Val(C_Script)

            ElseIf C_Script = "m_pre_pulse_sat" Then
                Call Advance_Script()
                If Val(C_Script) = 1 Then
                    pre_pulse_light = pre_pulse_light Or 1
                Else
                    pre_pulse_light = pre_pulse_light And (&HFFFF - 1)
                End If
                Jabber_jabber.List1.Items.Add("pre_pulse_light = " + Str(pre_pulse_light))


            ElseIf C_Script = "m_pre_pulse" Or C_Script = "pre_pulse" Then
                Dim mm As Integer
                For mm = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    ' Pre_Pulse(Current_Protocol, mm) = (C_Script)
                    'Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    If C_Script = "0" Then
                        Pre_Pulse(Current_Protocol, mm) = 0
                    Else
                        Pre_Pulse(Current_Protocol, mm) = 1
                        Pre_Pulse_Time(Current_Protocol, mm) = C_Script
                    End If

                Next mm

            ElseIf C_Script = "m_pre_pulse_delay" Or C_Script = "pre_delay" Then
                Dim mm As Integer
                For mm = 1 To M_Number_Loops(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Pre_Delay(Current_Protocol, mm) = C_Script
                Next mm

            ElseIf C_Script = "m_detector_gain" Then

                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Primary_Gain(Current_Protocol, M_Temp) = CShort(C_Script)
                Next M_Temp
                'Dim ix As Short
                'Jabber_jabber.gain_data_list.Items.Add("sample gains assigned from script m_detector_gain: ")
                'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '    Jabber_jabber.gain_data_list.Items.Add(M_Primary_Gain(Current_Protocol, ix))
                'Next

            ElseIf C_Script = "m_set_auto_gain" Then
                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    m_set_auto_gain(Current_Protocol, M_Temp) = CShort(C_Script)
                    Jabber_jabber.gain_data_list.Items.Add(m_set_auto_gain(Current_Protocol, M_Temp))
                Next M_Temp


                '
            ElseIf C_Script = "m_reference_gain" Then
                'Jabber_jabber.gain_data_list.Items.Add("ref gains during m_reference_gain: ")
                For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    Call Advance_Script()

                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    M_Reference_Gain(Current_Protocol, M_Temp) = CShort(C_Script)
                    Jabber_jabber.gain_data_list.Items.Add(M_Reference_Gain(Current_Protocol, M_Temp))
                Next M_Temp

                'Dim ix As Short
                'Jabber_jabber.gain_data_list.Items.Add("ref gains assigned from script m_reference_gain: ")
                'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '    Jabber_jabber.gain_data_list.Items.Add(M_Reference_Gain(Current_Protocol, ix))
                'Next


                ' ElseIf C_Script = "m_set_auto_gain" Then

                '     For M_Temp = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '         Call Advance_Script()
                '         'Script_Counter = Script_Counter + 1 'increment the script counter
                '         m_set_auto_gain(Current_Protocol, M_Temp) = CShort(C_Script)
                '         'Jabber_jabber.gain_data_list.Items.Add(m_set_auto_gain(Current_Protocol, M_Temp))
                '     Next M_Temp

                'NEW COMMANDS: 2010.09.07   9/7/2010
                ' runs the protocol saved in the given memory index
                ' usage: 
                '   current_protocol(2)
                '   upload_protocol
                '   ...
                '   run_protocol(2)
                '   ' the above loads protocol 2 into memory slot 2 of the FPGA
                '   ' and then runs protocol 2 from memory slot 2
            ElseIf C_Script = "run_protocol" Or C_Script = "run_protocol_from_memory" Or C_Script = "run" Then
                Dim memory As Integer
                Dim temporary_protocol As Short
                Call Advance_Script()
                memory = CInt(C_Script)
                temporary_protocol = Val(C_Script)
                If ForReal Then
                    Times_Averaged(Current_Trace) = Times_Averaged(Current_Trace) + 1
                    Start_Time = VB.Timer() - Zero_Time

                    Call Multi_Trace_from_FPGA(temporary_protocol, memory, M_Number_Loops, M_Intensity, _
                                               M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, _
                                               M_choke_actinic, L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, _
                                               Measuring_Pulse_Duration, Gain, In_Channel, _
                                               Ref_Channel_Number, _
                                               points, data_time, Xe_Flash, q_switch, M_Far_Red, M_Blue_Actinic, Pre_Pulse, _
                                               Pre_Pulse_Time, Pre_Delay, M_Take_Data)


                    Call Dirk_Save(data_time, points, File_Name(Current_Trace), Save_Mode(Current_Trace), _
                                   Times_Averaged(Current_Trace), Time_Mode(Current_Trace), Start_Time, _
                                   M_Number_Measuring_Lights(Current_Protocol), Baseline_Start(Current_Protocol), _
                                   Baseline_End(Current_Protocol), In_Channel(Current_Protocol), _
                                   Ref_Channel_Number(Current_Protocol), invert_raw, 1)
                    Me.Text = "SAVED:" & File_Name(Current_Trace)
                End If
            ElseIf C_Script = "upload_protocol" Or C_Script = "upload" Then
                ' programs the current protocol into the FPGA memory under the given index
                ' usage: 
                '   current_protocol(2)
                '   upload_protocol
                '   ' the above loads protocol 2 into memory slot 2 of the FPGA
                '   ' it is simplest to upload into the same memory slot as the protocol number

                Dim memory As Integer
                memory = CInt(Current_Protocol)
                If ForReal Then
                    Call SaveToChrisBlaster(Current_Protocol, memory, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, _
                                             M_Measuring_Light, M_choke_actinic, L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
                                              Gain, In_Channel, Ref_Channel_Number, points, data_time, Xe_Flash, q_switch, M_Far_Red, M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, M_Take_Data)

                    System.Windows.Forms.Application.DoEvents()
                End If

            ElseIf C_Script = "m_trace" Or C_Script = "m_trace_r" Then
                If ForReal Then


                    System.Windows.Forms.Application.DoEvents()




                    Times_Averaged(Current_Trace) = Times_Averaged(Current_Trace) + 1

                    if save_mode(current_trace) = Sequential then
                        file_name_sequence(Current_Trace) = file_name_sequence(Current_Trace) + 1
                        File_Name(Current_Trace) = Trim(Base_File_Name & VB6.Format(file_name_sequence(Current_Trace), "_0000_") & VB6.Format(Current_Trace, "0000") & ".dat")
                        'Aux_File_Name(i) = Trim(Base_File_Name & "_aux_" & VB6.Format(i, "0000") & ".dat")
                        Aux_File_Name(Current_Trace) = Trim(Base_File_Name & "_aux_" & VB6.Format(file_name_sequence(Current_Trace), "_0000_") & VB6.Format(Current_Trace, "0000") & ".dat")
                        Times_Averaged(Current_Trace) = 1
                    end if 



                    If Blocking_Filter(Current_Protocol) <> Last_Servo Then

                        If Blocking_Filter(Current_Protocol) = "blue" Then
                            SERVOset = Blue_Filter

                        ElseIf Blocking_Filter(Current_Protocol) = "ir" Then
                            SERVOset = IR_Filter

                        End If

                    End If
                    Last_Servo = Blocking_Filter(Current_Protocol)

                    Start_Time = VB.Timer() - Zero_Time 'save the start time of the trace

                    '                    Call Multi_Trace(Current_Protocol, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, M_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, Gain(Current_Protocol), In_Channel(Current_Protocol), points, data_time)

                    Call Multi_Trace(Current_Protocol, M_Number_Loops, M_Intensity, M_Number_Pulses, _
                                     M_Number_Measuring_Lights, M_Measuring_Light, M_choke_actinic, _
                                     L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
                                     Gain, In_Channel, Ref_Channel_Number, _
                                    points, data_time, Xe_Flash, _
                                     q_switch, M_Far_Red, M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, M_Take_Data)

                    If spec_mode = True Then
                        Dim irrr As Integer
                        For irrr = 0 To points
                            data_time(irrr) = start_lambda + irrr * delta_lambda
                            List1.Items.Add(data_time(irrr))
                        Next
                    End If


                    'List1.Items.Add("channels = " + Str(In_Channel(Current_Protocol)) + ", " + Str(Ref_Channel_Number(Current_Protocol)))
                    Call Dirk_Save(data_time, points, File_Name(Current_Trace), Save_Mode(Current_Trace), Times_Averaged(Current_Trace), Time_Mode(Current_Trace), Start_Time, M_Number_Measuring_Lights(Current_Protocol), Baseline_Start(Current_Protocol), Baseline_End(Current_Protocol), In_Channel(Current_Protocol), Ref_Channel_Number(Current_Protocol), invert_raw, 1)
                    'List1.Items.Add("aux channels = " + Str(Aux_In_Channel(Current_Protocol)) + ", " + Str(Aux_Ref_Channel(Current_Protocol)))


                    Call Dirk_Save(data_time, points, Aux_File_Name(Current_Trace), Save_Mode(Current_Trace), Times_Averaged(Current_Trace), Time_Mode(Current_Trace), Start_Time, M_Number_Measuring_Lights(Current_Protocol), Baseline_Start(Current_Protocol), Baseline_End(Current_Protocol), Aux_In_Channel(Current_Protocol), Aux_Ref_Channel(Current_Protocol), invert_raw, 0)

                    Me.Text = "SAVED:" & File_Name(Current_Trace)
                End If

            ElseIf C_Script = "offset" Then

                'Call servofrm.SETTHESERVO(1, (Wheel(Wheel_Position(Current_Protocol))))

                'Times_Averaged(Current_Trace) = Times_Averaged(Current_Trace) + 1

                If Blocking_Filter(Current_Protocol) <> Last_Servo Then

                    If Blocking_Filter(Current_Protocol) = "blue" Then
                        SERVOset = Blue_Filter
                        'Call servofrm.SETTHESERVO(0, SERVOset%)
                        'Call Hold_on_There(3)
                    ElseIf Blocking_Filter(Current_Protocol) = "ir" Then
                        SERVOset = IR_Filter
                        'Call servofrm.SETTHESERVO(0, SERVOset%)
                        'Call Hold_on_There(3)
                    End If

                End If
                Last_Servo = Blocking_Filter(Current_Protocol)
                If ForReal Then
                    Start_Time = VB.Timer() - Zero_Time 'save the start time of the trace
                    M_Number_Loops(0) = 1
                    M_Intensity(0, 1) = M_Intensity(Current_Protocol, 1)
                    M_Number_Pulses(0, 1) = 100
                    M_Number_Measuring_Lights(0) = 1
                    M_Number_Measuring_Lights(1) = 1
                    M_Measuring_Light(0, 1) = 0
                    M_choke_actinic(0, 1) = 0
                    M_Take_Data(0, 1) = 1
                    M_Measuring_Interval(0, 1) = "10m"
                    L_Measuring_Interval(0, 1, 1) = "10m"
                    M_Primary_Gain(0, 1) = M_Primary_Gain(Current_Protocol, 1)

                    Measuring_Pulse_Duration(0, 1) = Measuring_Pulse_Duration(Current_Protocol, 1)
                    Gain(0) = Gain(Current_Protocol)
                    In_Channel(0) = In_Channel(Current_Protocol)

                    '                    Call Multi_Trace(0, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, M_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, Gain(Current_Protocol), In_Channel(Current_Protocol), points, data_time)
                    Call Multi_Trace(0, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, _
                                     M_choke_actinic, L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
                                     Gain, In_Channel, Ref_Channel_Number, points, data_time, Xe_Flash, q_switch, _
                                     M_Far_Red, M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, M_Take_Data)
                    offset = 0

                    For i = 0 To points - 1
                        offset = offset + Data_Volts(i, In_Channel(Current_Protocol))
                    Next i

                    offset = offset / points
                    If ForReal Then List1.Items.Add("offset: " & offset)
                    'Scrolling text box fix: make the box jump to the bottom after adding an item
                    List1.TopIndex = List1.Items.Count - 1


                    '    Call Multi_Trace(0,
                    ', 1, M_Intensity(), 100, 1, 0, "10m", M_Primary_Gain(), Measuring_Pulse_Duration(), Gain(Current_Protocol), In_Channel(Current_Protocol), points&, data_time(), Data_Volts())


                    '    Call Dirk_Save(data_time(), Data_Volts(), points&, File_Name(Current_Trace), Save_Mode(Current_Trace), Times_Averaged(Current_Trace), Time_Mode(Current_Trace), Start_Time, M_Number_Measuring_Lights(Current_Protocol), Baseline_Start(Current_Protocol), Baseline_End(Current_Protocol))




                    'Dim M_Primary_Gain() As Integer
                End If
            ElseIf C_Script = "clear_gain_window" Then

                Jabber_jabber.gain_data_list.Items.Clear()


            ElseIf C_Script = "auto_gain" Then

                'Dim ix As Integer
                'Jabber_jabber.gain_data_list.Items.Add("before: ")
                'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '    Jabber_jabber.gain_data_list.Items.Add(M_Primary_Gain(Current_Protocol, ix))
                'Next
                'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                '    Jabber_jabber.gain_data_list.Items.Add(M_Reference_Gain(Current_Protocol, ix))
                'Next
                Dim temp_current_trace = Current_Trace


                If ForReal Then
                    'Dim auto_gain_output As Object
                    auto_gain_n(Current_Protocol, M_Number_Measuring_Lights, M_Number_Loops, M_Number_Pulses, _
                           Current_Intensity, M_Intensity, M_Take_Data, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
                            L_Measuring_Interval, Ref_Channel_Number, In_Channel, M_Measuring_Light, _
                            M_choke_actinic, Gain, points, data_time, Xe_Flash, q_switch, M_Far_Red, _
                            M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, File_Name, exclude_from_gain)

                    'Jabber_jabber.gain_data_list.Items.Add("after: ")
                    'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    '    Jabber_jabber.gain_data_list.Items.Add(M_Primary_Gain(Current_Protocol, ix))
                    'Next
                    'For ix = 1 To M_Number_Measuring_Lights(Current_Protocol)
                    '    Jabber_jabber.gain_data_list.Items.Add(M_Reference_Gain(Current_Protocol, ix))
                    'Next


                    'Dim t(,) As Short = M_Primary_Gain
                    'Dim tt(,) As Short = M_Reference_Gain
                    'Dim ttt(,) As String = Measuring_Pulse_Duration
                    Current_Trace = temp_current_trace
                End If

            ElseIf C_Script = "in_channel" Or C_Script = "meas_in_channel" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                In_Channel(Current_Protocol) = Val(C_Script)

                'finis
            ElseIf C_Script = "aux_in_channel" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Aux_In_Channel(Current_Protocol) = Val(C_Script)

                'finis
            ElseIf C_Script = "ref_in_channel" Or C_Script = "ref_channel_number" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Ref_Channel_Number(Current_Protocol) = Val(C_Script)

            ElseIf C_Script = "aux_ref_in_channel" Or C_Script = "aux_ref_channel" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Aux_Ref_Channel(Current_Protocol) = Val(C_Script)
            ElseIf C_Script = "spec_mode" Then ' set up the instument in spectro mode
                spec_mode = True

            ElseIf C_Script = "kin_mode" Then ' set up the instrument in kinetics mode
                spec_mode = False

            ElseIf C_Script = "start_lambda" Then ' set up the beginning wavelength for a spec trace
                Call Advance_Script()
                start_lambda = Val(C_Script)

            ElseIf C_Script = "delta_lambda" Then ' set up the number of nanometers to advance for each measuring pulse for a spec trace
                Call Advance_Script()
                delta_lambda = Val(C_Script)
                'finis

            ElseIf C_Script = "servo_positions" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                current_servo = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                number_servo_positions = Val(C_Script)
                For ii = 1 To number_servo_positions
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    servo_position(current_servo, ii) = Val(C_Script)
                Next ii


            ElseIf C_Script = "xe_intensity" Then


                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Xe_intensity_value = Val(C_Script)
                If ForReal Then List1.Items.Add("xe_intensity" & ":" & Xe_intensity_value)
                'Scrolling text box fix: make the box jump to the bottom after adding an item
                List1.TopIndex = List1.Items.Count - 1

                Call set_xe_intensity(Xe_intensity_value)
                ' ElseIf C_Script = "mobile_detector" Then
                '     Call Advance_Script()
                '     If C_Script = "fluor" Then
                '         tempi = 2
                '     ElseIf C_Script = "df" Then
                '         tempi = 0
                '     ElseIf C_Script = "absorbance" Then
                '         tempi = 1
                '     End If
                '     Call mobile_detector(tempi)
            ElseIf C_Script = "f_shutter" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then
                    'Call fluorescence_shutter(C_Script)
                    If new_firgelli = True Then


                        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)

                        For i = 1 To Val(C_Script)

                            Hold_on_There(0.001)

                            D_Out_Mask = (D_Out_Mask Or 64)
                            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
                            Label12.Text = i & "  " & D_Out_Mask

                            Hold_on_There(0.001)

                            D_Out_Mask = (D_Out_Mask And 191)
                            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)


                        Next
                    Else
                        Call fluorescence_shutter(C_Script)


                    End If
                    Call Hold_on_There(3)
                Else

                    If IsNumeric(C_Script) Then
                        ' input is fine
                    Else
                        'input is not fine
                        ShowScriptErrorMessage(Script_Counter, Script, "Value must be between 0 and 255")
                    End If
                End If
            ElseIf C_Script = "stir" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call stir_control(C_Script)

            ElseIf C_Script = "subtract" Then

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace1 = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace2 = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace3 = Val(C_Script)
                If ForReal Then Call subtract_traces(File_Name(trace1), File_Name(trace2), File_Name(trace3))

            ElseIf C_Script = "dice_average" Then

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace1 = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace2 = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                dice_length = Val(C_Script)
                If ForReal Then Call dice_average(File_Name(trace1), File_Name(trace2), dice_length)

            ElseIf C_Script = "delta_a" Then
                If ForReal Then
                    'If use_ref_channel = False Then
                    Call Delta_A_calculate_smooth_reference(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)), 0, alternate_reference_file_name(Current_Trace))

                    'Call Delta_A_calculate(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)))
                    'Else
                    'Call Delta_A_WITH_ref(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)))
                    'End If
                End If
            ElseIf C_Script = "delta_a_s" Then
                If ForReal Then
                    Call Advance_Script()
                    Dim smooth_window As Short = Val(C_Script)
                    Call Delta_A_calculate_smooth_reference(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)), smooth_window, alternate_reference_file_name(Current_Trace))
                End If

            ElseIf C_Script = "ratio_calc" Then
                If ForReal Then
                    'If use_ref_channel = False Then
                    Call Ratio_Calculate(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)))
                    'Else
                    'Call Delta_A_WITH_ref(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)))
                    'End If
                End If
            ElseIf C_Script = "normalize_to_integration_time" Then
                If ForReal Then
                    'If use_ref_channel = False Then
                    Call normalize_to_integration_time(File_Name(Current_Trace), Measuring_Pulse_Duration, Current_Protocol)
                    'Else
                    'Call Delta_A_WITH_ref(File_Name(Current_Trace), (Baseline_Start(Current_Protocol)), (Baseline_End(Current_Protocol)))
                    'End If
                End If
            ElseIf C_Script = "zero_point" Then
                Call Advance_Script()
                If ForReal Then Call Zero_Point(File_Name(Current_Trace), Val(C_Script))

            ElseIf C_Script = "extract_points_delta" Then
                If ForReal Then
                    Dim extract_index, extract_file_number, extract_this_point, extract_from_this_number_of_traces As Integer
                    Dim extracted_value As Single
                    FileOpen(6, File_Name(Current_Trace), OpenMode.Output, OpenAccess.Default, OpenShare.Shared) 'open current trace file for output
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    extract_this_point = Val(C_Script)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    extract_from_this_number_of_traces = Val(C_Script)
                    For extract_index = 1 To extract_from_this_number_of_traces
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1
                        extracted_value = extract_point(File_Name(C_Script), extract_this_point)
                        extract_file_number = Val(C_Script)
                        If ForReal Then List1.Items.Add("extracted: " & extracted_value)
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        PrintLine(6, extracted_value)

                    Next extract_index
                    FileClose(6)
                    'Threading.Thread.Sleep(10)

                End If
            ElseIf C_Script = "extract_points_raw" Then
                If ForReal Then
                    Dim extract_index, extract_file_number, extract_this_point, extract_from_this_number_of_traces As Integer
                    Dim extracted_value As Single
                    FileOpen(6, File_Name(Current_Trace), OpenMode.Output, OpenAccess.Default, OpenShare.Shared) 'open current trace file for output
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    extract_this_point = Val(C_Script)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    extract_from_this_number_of_traces = Val(C_Script)
                    For extract_index = 1 To extract_from_this_number_of_traces
                        Call Advance_Script()
                        'Script_Counter = Script_Counter + 1
                        extracted_value = extract_point_raw(File_Name(C_Script), extract_this_point)
                        extract_file_number = Val(C_Script)
                        If ForReal Then List1.Items.Add("extracted: " & extracted_value)
                        'Scrolling text box fix: make the box jump to the bottom after adding an item
                        List1.TopIndex = List1.Items.Count - 1
                        PrintLine(6, extracted_value)

                    Next extract_index
                    FileClose(6)
                    'Threading.Thread.Sleep(10)

                End If
            ElseIf C_Script = "load_file" Then  'load a trace into memory

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                File_Name(Current_Trace) = C_Script


            ElseIf C_Script = "set_servo" Then

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                current_servo = Val(C_Script)

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                current_servo_position = Val(C_Script)

                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call set_servo(current_servo, servo_position(current_servo, current_servo_position))

            ElseIf C_Script = "simmer" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                simmer_number = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                q_switch_time = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                simmer_q_switch = Val(C_Script)
                If ForReal Then Call Simmer(simmer_number, (Str(q_switch_time) & "u"), simmer_q_switch)

            ElseIf C_Script = "number_traces" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Number_Traces = Val(C_Script) '

                ReDim Times_Averaged(Number_Traces)
                ReDim Time_Mode(Number_Traces)
                ReDim trace_note(Number_Traces)
                ReDim trace_label(Number_Traces)

                ReDim smooth(Number_Traces)

                ReDim Save_Mode(Number_Traces)
                ReDim file_name_sequence(Number_Traces)

                ReDim File_Name(Number_Traces)
                ReDim alternate_reference_file_name(Number_Traces)
                ReDim Aux_File_Name(Number_Traces)
                ReDim Trace_protocol(Number_Traces)
                '    m Baseline_Start(Number_Traces)
                '   ReDim Baseline_End(Number_Traces)

                For i = 0 To Number_Traces ' put in some default values
                    Save_Mode(i) = Average_Into
                    file_name_sequence(i)=0
                    Times_Averaged(i) = 0
                    Trace_protocol(i) = -1
                    File_Name(i) = Trim(Base_File_Name & VB6.Format(i, "0000") & ".dat")
                    Aux_File_Name(i) = Trim(Base_File_Name & "_aux_" & VB6.Format(i, "0000") & ".dat")
                    alternate_reference_file_name(i) = ""
                    'List1.AddItem File_Name(i)
                    Time_Mode(i) = "from_zero"

                Next i



            ElseIf C_Script = "purge_file" Then
                If ForReal Then
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    'Recovery_Measuring_Interval(Current_Protocol) = c_script
                    FileOpen(6, File_Name(Current_Trace), OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
                    PrintLine(6, "")
                    FileClose(6)
                    ' Threading.Thread.Sleep(10)

                End If
            ElseIf C_Script = "current_trace" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                If C_Script = "v" Then
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Array_Number = Val(C_Script)
                    'UPGRADE_WARNING: Couldn't resolve default property of object Array_Value(). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                    Current_Trace = Val(Array_Value(Array_Number, Array_Index(Array_Number)))
                Else
                    Current_Trace = Val(C_Script) 'set wait time
                End If
                If Trace_protocol(Current_Trace) >= 0 Then  'if we have set a trace_protocol, then set current_protocol to it
                    Current_Protocol = Trace_protocol(Current_Trace)
                End If
            ElseIf C_Script = "trace_protocol" Then
                Call Advance_Script()
                Trace_protocol(Current_Trace) = Val(C_Script) 'set wait time


            ElseIf C_Script = "time_mode" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Time_Mode(Current_Trace) = C_Script 'set wait time
                ' time_mode="from_zero", or ""sequential"
            ElseIf C_Script = "trace_note" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace_note(Current_Trace) = C_Script '
                ' time_mode="from_zero", or ""sequential"

            ElseIf C_Script = "trace_label" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                trace_label(Current_Trace) = C_Script '
                ' time_mode="from_zero", or ""sequential"

            ElseIf C_Script = "number_protocols" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Number_Protocols = Val(C_Script) 'set wait time



                ReDim M_Number_Loops(Number_Protocols)


                ReDim M_Intensity(Number_Protocols, Max_Number_Loops)
                ReDim M_Far_Red(Number_Protocols, Max_Number_Loops)
                ReDim M_Blue_Actinic(Number_Protocols, Max_Number_Loops)
                ReDim M_Number_Pulses(Number_Protocols, Max_Number_Loops)
                ReDim M_Number_Measuring_Lights(Number_Protocols)
                ReDim Pre_Pulse(Number_Protocols, Max_Number_Loops)
                ReDim Pre_Delay(Number_Protocols, Max_Number_Loops)
                ReDim Pre_Pulse_Time(Number_Protocols, Max_Number_Loops)

                'Measuring_Interval(Current_Protocol, lightindex + 1, loopindex + 1)
                'Measuring_Light(Current_Protocol, lightindex + 1)
                'choke_actinic(Current_Protocol, lightindex + 1)
                'Primary_Gain(Current_Protocol, lightindex + 1)
                'reference_gain(Current_Protocol, lightindex + 1)

                ReDim M_Take_Data(Number_Protocols, Max_Number_Loops)
                ReDim M_Measuring_Light(Number_Protocols, Max_Number_Lights)
                ReDim M_Measuring_Light_Names(Number_Protocols, Max_Number_Lights)

                '                ReDim M_Wasp(Number_Protocols, Max_Number_Lights)

                ReDim m_plot_specific_traces(Number_Protocols + 1, Max_Number_Lights + 1)

                ReDim M_choke_actinic(Number_Protocols, Max_Number_Lights)

                ReDim M_Measuring_Interval(Number_Protocols, Max_Number_Lights)
                ReDim M_Reference_Gain(Number_Protocols, Max_Number_Lights)

                ReDim L_Measuring_Interval(Number_Protocols, Max_Number_Lights, Max_Number_Loops)
                ReDim M_Primary_Gain(Number_Protocols, Max_Number_Lights)
                ReDim m_set_auto_gain(Number_Protocols, Max_Number_Lights)

                ReDim Measuring_Pulse_Duration(Number_Protocols, Max_Number_Lights)

                ReDim Xe_Flash(Number_Protocols, Max_Number_Loops)
                ReDim q_switch(Number_Protocols, Max_Number_Loops)
                ReDim Baseline_Start(Number_Protocols)
                ReDim Baseline_End(Number_Protocols)
                ReDim Detector_Gain(Number_Protocols)
                ReDim In_Channel(Number_Protocols)
                ReDim Aux_In_Channel(Number_Protocols)
                ReDim Ref_Channel_Number(Number_Protocols)
                ReDim Aux_Ref_Channel(Number_Protocols)

                ReDim Wavelength(Number_Protocols)
                ReDim Gain(Number_Protocols)
                ReDim Blocking_Filter(Number_Protocols)
                ReDim Wheel_Position(Number_Protocols)

                'ReDim Dirk_Measuring_Interval(Number_Protocols)
                'ReDim Number_Dirk_Repeats(Number_Protocols)
                ReDim Baseline_Measuring_Interval(Number_Protocols)
                ReDim Recovery_Measuring_Interval(Number_Protocols)
                ReDim Baseline_Number_Points(Number_Protocols)
                ReDim Recovery_Number_Points(Number_Protocols)
                ReDim Dirk_Number_Points(Number_Protocols)
                Dim Actinic_Intensity(Number_Protocols) As Object
                'ReDim Baseline_Actinic_Intensity(Number_Protocols)
                'ReDim Recovery_Actinic_Intensity(Number_Protocols)
                'ReDim Dirk_Actinic_Intensity(Number_Protocols)
                'ReDim Measuring_Mode(Number_Protocols)
                ReDim Give_Flash(Number_Protocols)

                'set some default values


                For i = 0 To Number_Protocols
                    For ii = 0 To Max_Number_Loops
                        M_Take_Data(i, ii) = 1

                    Next
                Next
                For i = 0 To Number_Protocols
                    For ii = 0 To Max_Number_Lights
                        M_choke_actinic(i, ii) = 0
                        m_set_auto_gain(i, ii) = 1 'the default value is 1, i.e. set it
                        m_plot_specific_traces(i, ii) = 1
                    Next
                Next


                For i = 1 To Number_Protocols
                    M_Number_Measuring_Lights(i) = 1
                    Detector_Gain(i) = 4
                    Wavelength(i) = "520"
                    Gain(i) = BIP10VOLTS
                    Blocking_Filter(i) = "blue"
                    Wheel_Position(i) = 2
                    In_Channel(i) = 0
                    Ref_Channel_Number(i) = 1
                    Aux_In_Channel(i) = 2
                    Aux_Ref_Channel(i) = 3


                    For nlight As Integer = 1 To Max_Number_Lights
                        Measuring_Pulse_Duration(i, nlight) = "20u"
                        M_Measuring_Light_Names(i, nlight - 1) = Str(nlight)
                        'M_Wasp(i, nlight) = 0
                    Next nlight

                    'UPGRADE_WARNING: Couldn't resolve default property of object Actinic_Intensity(i). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                    Actinic_Intensity(i) = 7
                    Give_Flash(i) = False
                Next i

                Current_Protocol = 1

            ElseIf C_Script = "current_protocol" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Current_Protocol = Val(C_Script) 'set wait time

            ElseIf C_Script = "blocking_filter" Then
                Call Advance_Script()
                'Blocking_Filter(Current_Protocol) = C_Script "blocking_filter" used to be set as a servo property
                ' Now it is the same as f_shutter
                If ForReal Then
                    Call fluorescence_shutter(C_Script)
                Else
                    If C_Script = "open" Or C_Script = "1" Or C_Script = "close" Or C_Script = "0" Then
                        ' input is fine
                    Else
                        'input is not fine
                        ShowScriptErrorMessage(Script_Counter, Script, "Value must be ""1"", ""open"", ""0"", or ""close"".")
                    End If
                End If
            ElseIf C_Script = "home_wheel" Then
                If ForReal Then Call home_wheel()
            ElseIf C_Script = "wheel_position" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                temp_s = Val(C_Script)
                If ForReal Then Call set_wheel(temp_s)
            ElseIf C_Script = "give_flash" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Give_Flash(Current_Protocol) = Val(C_Script)
            ElseIf C_Script = "nmp" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                nmp = Val(C_Script) 'set wait time

                'ElseIf C_Script = "measuring_mode" Then
                '        Call Advance_Script()
                '        'Script_Counter = Script_Counter + 1 'increment the script counter
                '        Measuring_Mode(Current_Protocol) = C_Script 'set wait time
                '        If Measuring_Mode(Current_Protocol) = "520" Then
                '            Wavelength(Current_Protocol) = "520"
                '            Blocking_Filter(Current_Protocol) = "blue"
                '        ElseIf Measuring_Mode(Current_Protocol) = "l2" Then
                '            Wavelength(Current_Protocol) = "???"
                '            Blocking_Filter(Current_Protocol) = "blue"

                '        ElseIf Measuring_Mode(Current_Protocol) = "820" Then
                '            Wavelength(Current_Protocol) = "820"
                '            Blocking_Filter(Current_Protocol) = "ir"
                '        ElseIf Measuring_Mode(Current_Protocol) = "fluor" Then
                '            Wavelength(Current_Protocol) = "520"
                '            Blocking_Filter(Current_Protocol) = "ir"
                '        ElseIf Measuring_Mode(Current_Protocol) = "sat" Then
                '            Wavelength(Current_Protocol) = "520"
                '            Blocking_Filter(Current_Protocol) = "ir"
                '        ElseIf Measuring_Mode(Current_Protocol) = "satw" Then
                '            Wavelength(Current_Protocol) = "520"
                '            Blocking_Filter(Current_Protocol) = "ir"
                '        ElseIf Measuring_Mode(Current_Protocol) = "blue_artifact" Then
                '            Wavelength(Current_Protocol) = "none"
                '            Blocking_Filter(Current_Protocol) = "blue"
                '        ElseIf Measuring_Mode(Current_Protocol) = "actinic_intensity" Then
                '            Wavelength(Current_Protocol) = "none"
                '            Blocking_Filter(Current_Protocol) = "blue"
                '            In_Channel(Current_Protocol) = 2
                '        End If

            ElseIf C_Script = "detector_gain" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Detector_Gain(Current_Protocol) = CShort(C_Script) 'set wait time

            ElseIf C_Script = "gain" Or C_Script = "adc_gain" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                temp = C_Script 'set wait time
                If temp = "bip10volts" Then
                    Gain(Current_Protocol) = BIP10VOLTS
                ElseIf temp = "bip5volts" Then
                    Gain(Current_Protocol) = BIP5VOLTS
                ElseIf temp = "bip2volts" Then
                    Gain(Current_Protocol) = BIP2VOLTS
                ElseIf temp = "bip1volts" Then
                    Gain(Current_Protocol) = BIP1VOLTS
                    'ElseIf temp = "bippt25volts" Then
                    '   Gain(Current_Protocol) = BIPPT25VOLTS

                Else
                    Gain(Current_Protocol) = BIP10VOLTS 'if something screws up, give +/5 5v
                End If


            ElseIf C_Script = "message" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                temp_string = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                If C_Script = "v" Then
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Array_Number = Val(C_Script)

                    temp_string = temp_string & Str(Array_Value(Array_Number, Array_Index(Array_Number)))
                Else
                    temp_string = temp_string & C_Script
                End If
                'Me.Text = temp_string

            ElseIf C_Script = "mass_flow_2" Or C_Script = "mass_flow_room" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                Intensity_Form.TextBox4.Text = Val(C_Script)
                Call Intensity_Form.mass_flow_set()

            ElseIf C_Script = "mass_flow_3" Or C_Script = "mass_flow_scrub" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                Intensity_Form.TextBox5.Text = Val(C_Script)
                Call Intensity_Form.mass_flow_set()

            ElseIf C_Script = "intensity" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                If C_Script = "v" Then
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    Array_Number = Val(C_Script)

                    Current_Intensity = Array_Value(Array_Number, Array_Index(Array_Number))
                Else
                    Current_Intensity = Val(C_Script) 'set current intensity
                End If
                If ForReal Then
                    Call Out_Intensity(Current_Intensity, Blue_Actinic, Far_Red, 0)

                    ProgressBar3.Value = Current_Intensity
                    actinic_label.Text = Current_Intensity
                End If
            ElseIf C_Script = "blue_actinic" Or C_Script = "flash" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Blue_Actinic = Val(C_Script) 'set current intensity
                If ForReal Then Call Out_Intensity(Current_Intensity, Blue_Actinic, Far_Red, 0)

            ElseIf C_Script = "leaf_actinic" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Blue_Actinic = Val(C_Script) 'set current intensity
                If ForReal Then Call Out_Intensity(Current_Intensity, Blue_Actinic, Far_Red, 0)
            ElseIf C_Script = "far_red" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Far_Red = Val(C_Script) 'set current intensity
                If ForReal Then Call Out_Intensity(Current_Intensity, Blue_Actinic, Far_Red, 0)

            ElseIf C_Script = "recovery_number_points" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Recovery_Number_Points(Current_Protocol) = Val(C_Script) 'set wait time

                'ElseIf C_Script = "baseline_actinic_intensity" Then
                '        Call Advance_Script()
                '        'Script_Counter = Script_Counter + 1 'increment the script counter

                '        If C_Script = "v" Then
                '            Call Advance_Script()
                '            'Script_Counter = Script_Counter + 1 'increment the script counter
                '            Array_Number = Val(C_Script)

                '            Baseline_Actinic_Intensity(Current_Protocol) = Array_Value(Array_Number, Array_Index(Array_Number))
                '        Else
                '            Baseline_Actinic_Intensity(Current_Protocol) = Val(C_Script) 'set wait time
                '        End If
                'ElseIf C_Script = "recovery_actinic_intensity" Then
                '        Call Advance_Script()
                '        'Script_Counter = Script_Counter + 1 'increment the script counter
                '        If C_Script = "v" Then
                '            Call Advance_Script()
                '            'Script_Counter = Script_Counter + 1 'increment the script counter
                '            Array_Number = Val(C_Script)
                '            'UPGRADE_WARNING: Couldn't resolve default property of object Array_Value(Array_Number, Array_Index()). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                '            Recovery_Actinic_Intensity(Current_Protocol) = Array_Value(Array_Number, Array_Index(Array_Number))
                '        Else
                '            Recovery_Actinic_Intensity(Current_Protocol) = Val(C_Script)
                '        End If
                'ElseIf C_Script = "dirk_actinic_intensity" Then
                '        Call Advance_Script()
                '        'Script_Counter = Script_Counter + 1 'increment the script counter
                '        If C_Script = "v" Then
                '            Call Advance_Script()
                '            'Script_Counter = Script_Counter + 1 'increment the script counter
                '            Array_Number = Val(C_Script)
                '            'UPGRADE_WARNING: Couldn't resolve default property of object Array_Value(Array_Number, Array_Index()). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                '            Dirk_Actinic_Intensity(Current_Protocol) = Array_Value(Array_Number, Array_Index(Array_Number))
                '        Else
                '            Dirk_Actinic_Intensity(Current_Protocol) = Val(C_Script) 'set wait time
                '        End If
            ElseIf C_Script = "no_plot" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                No_Plot_Data = Val(C_Script)

                'ElseIf C_Script = "smooth" Then
                ' not implemented

            ElseIf C_Script = "plot_file" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script m_plot_specific_traces
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (1), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "plot_clear" Then
                For i = 1 To 6
                    'List1.Items.Add(Graph_Name(i))
                    If ForReal Then Call Plot_File(rootdata_directory & "test.dat", (i), (i And 1), False, 0, Current_Protocol, m_plot_specific_traces)
                Next i

            ElseIf C_Script = "plot_raw" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (0), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "plot_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (1), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "plot_delta_linear" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (1), False, 1, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "plot_ref" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (2), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "add_plot_raw" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (0), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "add_plot_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (1), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "add_plot_delta_linear" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (1), True, 1, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "add_plot_ref" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (2), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)


            ElseIf C_Script = "calc_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Delta(File_Name(Current_Trace), File_Name(Val(C_Script)))
            ElseIf C_Script = "combine_json" Then

                combineJSON(Base_File_Name)

            ElseIf C_Script = "json" Then
                'Call Advance_Script()

                If ForReal Then Call saveJSON(File_Name(Current_Trace), Aux_File_Name(Current_Trace), Current_Protocol, _
                        M_Number_Loops, M_Intensity, _
                        M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, M_choke_actinic, L_Measuring_Interval, _
                            M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
                            Gain, In_Channel, Ref_Channel_Number, _
                            Xe_Flash, q_switch, M_Far_Red, M_Blue_Actinic, trace_note, trace_label, M_Measuring_Light_Names, Time_Mode, Start_Time, _
                            Current_Trace, protocol_label, current_lambda)

                'Pre_Pulse_Time

                ')

                '                saveJSON(ByRef I_File_Name As String, ByRef I0_File_Name As String, ByRef Current_Protocol As Short, ByRef Number_Loops() As Short, ByRef Intensity(,) _
                '                            As Short, ByRef Number_Pulses(,) As Short, ByRef Number_Measuring_Lights() As Short, _
                '                            ByRef Measuring_Light(,) As Short, ByRef choke_actinic(,) As Short, ByRef Measuring_Interval(,,) As String, _
                '                            ByRef Primary_Gain(,) As Short, ByRef reference_gain(,) As Short, ByRef Measuring_Pulse_Duration(,) As String, _
                '                            ByRef Gain() As Integer, ByRef In_Channel() As Short, ByRef Ref_channel_number() As Short, ByRef points As Integer, ByRef data_time() As Single, _
                '                            ByRef Xe_Flash(,) As Short, ByRef q_switch(,) As Short, ByVal M_Far_Red(,) As Short, ByVal M_Blue_actinic(,) _
                '                            As Short, ByVal pre_pulse(,) As Short, ByVal pre_pulse_time(,) As String, ByVal Pre_Delay(,) As String, _
                '                            ByVal m_take_data(,) As Short,
                '                            ByVal M_Wasp(,) As Short)


            ElseIf C_Script = "aux_plot_raw" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then
                    Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (0), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

                End If
                'If ForReal Then Call Plot_File(File_Name(Current_Trace), (Val(C_Script)), (0), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "aux_plot_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (1), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "aux_plot_delta_linear" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (1), False, 1, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "aux_plot_ref" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (2), False, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "aux_add_plot_raw" Or C_Script = "add_aux_plot_raw" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (0), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "aux_add_plot_delta" Or C_Script = "add_aux_plot_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (1), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)
            ElseIf C_Script = "aux_add_plot_delta_linear" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (1), True, 1, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "aux_add_plot_ref" Or C_Script = "add_aux_plot_ref" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Plot_File(Aux_File_Name(Current_Trace), (Val(C_Script)), (2), True, 0, Current_Protocol, m_plot_specific_traces) ', Graph_Name)

            ElseIf C_Script = "aux_calc_delta" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                If ForReal Then Call Delta(Aux_File_Name(Current_Trace), File_Name(Val(C_Script)))









            ElseIf C_Script = "sub_var" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                i = (sub_var(C_Script, 1))

            ElseIf C_Script = "save_mode" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                temp = C_Script
                If temp = "file_replace" Or temp = "overwrite" Then
                    Save_Mode(Current_Trace) = File_Replace
                ElseIf temp = "average_into" Then
                    Save_Mode(Current_Trace) = Average_Into
                ElseIf temp = "file_append" Then
                    Save_Mode(Current_Trace) = File_Append
                ElseIf temp = "sequential" Then
                    Save_Mode(Current_Trace) = Sequential
                End If

                'Global Const File_Replace = 0
                'Global Const Average_Into = 1
                'Global Const File_Append = 2
            ElseIf C_Script = "graph_labels" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter


                Graph_Name(1) = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Graph_Name(2) = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Graph_Name(3) = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Graph_Name(4) = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Graph_Name(5) = C_Script
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Graph_Name(6) = C_Script


            ElseIf C_Script = "array_index" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Array_Number = Val(C_Script)
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Array_Index(Array_Number) = Val(C_Script)

            ElseIf C_Script = "inc_array" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Array_Number = Val(C_Script) 'set wait time
                Array_Index(Array_Number) = Array_Index(Array_Number) + 1

            ElseIf C_Script = "dec_array" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter
                Array_Number = Val(C_Script) 'set wait time
                Array_Index(Array_Number) = Array_Index(Array_Number) - 1
            ElseIf VB.Right(C_Script, 1) = "|" Then
                ' indicates the start of a subroutine, ignore it
            ElseIf C_Script = "record_files" Then  'adds text to the note file
                Call Advance_Script()
                If Val(C_Script) > 0 Then
                    record_files = True
                Else
                    record_files = False
                End If

            ElseIf C_Script = "v_array" Then
                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                Array_Number = Val(C_Script)

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                'UPGRADE_WARNING: Couldn't resolve default property of object Array_Max_Index(Array_Number). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                Array_Max_Index(Array_Number) = Val(C_Script)

                'UPGRADE_WARNING: Couldn't resolve default property of object Array_Max_Index(Array_Number). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                For Array_Index(Array_Number) = 1 To Array_Max_Index(Array_Number)
                    Call Advance_Script()
                    'Script_Counter = Script_Counter + 1 'increment the script counter
                    'UPGRADE_WARNING: Couldn't resolve default property of object Array_Value(Array_Number, Array_Index()). Click for more: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
                    Array_Value(Array_Number, Array_Index(Array_Number)) = Val(C_Script)
                Next



            ElseIf C_Script = "att_sat" Then  ' change position of the fv/fv saturation pulse attenuation disk

                Call Advance_Script()
                'Script_Counter = Script_Counter + 1 'increment the script counter

                If Val(C_Script) = 1 Then
                    wt = 185
                Else
                    wt = 35
                End If

                'Call servofrm.SETTHESERVO(2, wt)

                'Else
                '    System.Windows.Forms.MessageBox.Show("Error:Command " & C_Script & " not recognized")
                '    Exit Sub
            Else ' script function not recognized

                Dim keepgoing As Boolean
                keepgoing = ShowScriptErrorMessage(Script_Counter, Script, "This command was not recognized.")
                If keepgoing = False Then Exit Sub
            End If

        End While
        If ForReal Then

            'If record_files = True Then  'adds text to the note file")
            '    note_text = note_text & "*******************************************************************" & ControlChars.CrLf
            '    note_text = note_text & "ALL FILE FILE NAMES IN SCRIPT" & ControlChars.CrLf


            '    For i = 0 To Number_Traces
            '        If Times_Averaged(i) > 0 Then
            '            note_text = note_text & File_Name(i) & ControlChars.CrLf & " averaged " & Times_Averaged(i) & " times. " & ControlChars.CrLf
            '            note_text = note_text & trace_note(i) & ControlChars.CrLf
            '        End If
            '    Next i
            '    note_text = note_text & "*******************************************************************" & ControlChars.CrLf
            '    FileOpen(7, Note_File_Name, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
            '    Write(7, note_text)
            '    Write(7, '\n')
            '    FileClose(7)
            '    'Threading.Thread.sleep(10)

            'End If
            If Halt_Script = True Then
                ulstat = cbStopBackground(BoardNum, AIFUNCTION)
                List1.Items.Add("HALT!")
                'Scrolling text box fix: make the box jump to the bottom after adding an item
                List1.TopIndex = List1.Items.Count - 1
            End If
spook:      Exit Sub

err_it:
            If Err.Number = 32755 Then
                Err.Clear()
                CommonDialog1Open.Reset()


                Resume spook
            Else

                ShowScriptFatalErrorMessage(Script_Counter, Script, "Error occured while reading script file: " & vbNewLine & Err.Number & ": " & Err.Description)
                'System.Windows.Forms.MessageBox.Show(Err.Number & " " & Err.Description & " Error in the script file around " & Script_Counter & " text = " & C_Script)
                Err.Clear()
                Close()

                Script_Counter = 0
                Resume spook
            End If
        Else ' end of test run
            List1.Items.Add("reached the end of the script file")
        End If

    End Sub


    Public Sub set_base_file_name()
        'On Error GoTo err_it3
        Dim prev As String = ""
        CommonDialog1Open.Title = "enter file name (& location) for base file name"
        CommonDialog1Save.Title = "enter file name (& location) for base file name"
        CommonDialog1Open.FileName = "*.dat"
        CommonDialog1Save.FileName = "*.dat"
        CommonDialog1Open.Title = "base file name to save"
        CommonDialog1Save.Title = "base file name to save"
        CommonDialog1Save.ShowDialog()
        CommonDialog1Open.FileName = CommonDialog1Save.FileName


        If CommonDialog1Open.FileName = "" Then
            Base_File_Name = "temp"
            Base_Base_File_Name = "temp"
        Else
            Base_Base_File_Name = CommonDialog1Open.FileName
            Base_File_Name = CommonDialog1Open.FileName
            List1.Items.Add("& " & Base_Base_File_Name)
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1

        End If
        Exit Sub
err_it3:

        Err.Clear()
        script_file = prev
        Script_Label.Text = script_file
        Resume get_out
get_out:

    End Sub
    Public Sub set_auto_file_name()
        'On Error GoTo err_it3
        Dim prev As String = ""
        CommonDialog1Open.Title = "enter file name (& location) for AUTORUN"
        CommonDialog1Save.Title = "enter file name (& location) for AUTORUN"
        CommonDialog1Open.FileName = "*.txt"
        CommonDialog1Save.FileName = "*.txt"
        CommonDialog1Open.Title = "AUTORUN.TXT in this folder"
        CommonDialog1Save.Title = "AUTORUN.TXT in this folder"
        CommonDialog1Save.ShowDialog()
        CommonDialog1Open.FileName = CommonDialog1Save.FileName


        If CommonDialog1Open.FileName = "" Then
            Auto_File_Name = "temp"

        Else
            Auto_File_Name = CommonDialog1Open.FileName
            'Base_File_Name = CommonDialog1Open.FileName
            List1.Items.Add("& " & Auto_File_Name)
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1

        End If
        Exit Sub
err_it3:

        Err.Clear()
        script_file = prev
        Script_Label.Text = script_file
        Resume get_out
get_out:

    End Sub
    Public Sub Advance_Script()
        Script_Counter = Script_Counter + 1 'advance the script counter
        'List1.Items.Clear()



        If VB.Left(Script(Script_Counter), 1) = "@" Then 'check to see if the next script line designates that we are inserting variable
            Dim j As Single
            Dim i, k As Integer
            Dim var As String

            var = VB.Right(Script(Script_Counter), (Len(Script(Script_Counter)) - 1))  'which variable are we inserting?
            'Label7.Text = var
            j = 0
            For i = 1 To number_variables
                If var = user_variable_name(i) Then
                    j = i
                End If
            Next i
            If j = 0 Then
                System.Windows.Forms.MessageBox.Show("Error:I cannot find that variable", "ERROR", MessageBoxButtons.OK, MessageBoxIcon.Error)
                Exit Sub
            Else
                Script_Counter = Script_Counter + 1
                If Script(Script_Counter) = "inc" Then
                    user_variable_index(j) = user_variable_index(i) + 1
                    k = user_variable_index(j)
                ElseIf Script(Script_Counter) = "dec" Then
                    user_variable_index(j) = user_variable_index(i) + 1
                    k = user_variable_index(j)
                Else
                    k = Val(Script(Script_Counter))
                    user_variable_index(j) = k
                End If
                'Script_Counter = Script_Counter + 1

                C_Script = (user_variable(j, k)) 'insert the value of the variable into the script
            End If
        Else
            C_Script = Script(Script_Counter)
        End If


    End Sub
    Public Function sub_var(ByVal var As String, ByVal ii As Integer) As Single
        Dim j As Single
        Dim i As Integer
        If VB.Left(var, 1) = "@" Then
            For i = 1 To number_variables
                If var = user_variable_name(i) Then
                    j = user_variable(i, ii)
                End If
            Next i
        Else
            j = Val(var)
        End If
        List1.Items.Add(var & " " & j)
        'Scrolling text box fix: make the box jump to the bottom after adding an item
        List1.TopIndex = List1.Items.Count - 1
        Return (j)
    End Function

    Private Sub Script_Label_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles Script_Label.Click
        On Error GoTo err_it2
        Dim prev As String
        prev = script_file
        'Dim script_file As String
        If prev <> "" Then
            CommonDialog1Open.FileName = prev
        Else
            CommonDialog1Open.FileName = "*.txt"
        End If
        'CommonDialog1Open.Filter = "Script Files|*.scr|Text Files|*.txt"
        CommonDialog1Open.Filter = "Text Files|*.txt"
        CommonDialog1Save.FileName = "*.txt"
        CommonDialog1Open.ShowDialog()
        '        CommonDialog1Save.FileName = CommonDialog1Open.FileName


        script_file = CommonDialog1Open.FileName
        Script_Label.Text = script_file
        Exit Sub
err_it2:
        Err.Clear()
        script_file = prev
        Script_Label.Text = script_file
        Resume get_out
get_out:
    End Sub

    Public Sub servoform_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs)
        'servofrm.Show()

    End Sub

    Public Sub Set_Actinic_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles Set_Actinic.Click
        Intensity_Form.Show()
        Intensity_Form.BringToFront()

        'Intensity_done = False

        'While Intensity_done = False
        ' System.Windows.Forms.Application.DoEvents()
        ' End While

        'Call Out_Intensity(Intensity_Pass, blue_actinic, far_red)

    End Sub

    Public Sub skip_wait_Click(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles skip_wait.Click
        Bug_Out = True

    End Sub

    Public Sub set_wheel(ByRef wp As Short)
        'Call servofrm.SETTHESERVO(1, wp)
    End Sub
    Private Sub ClearWindow(ByVal ThisGraph As Short)
        'this sub chooses the first picturebox as the default image for clearing
        'it changes all the pixels to black in the chosen box


    End Sub



    Private Sub Button1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Dim i As Integer
        'Dim x As Double
        Dim arrData(30, 0) As Single

        'Dim i, j As Integer
        'Dim rowLabelCount, columnLabelCount, rowCount As Integer
        'Dim columnCount, labelIndex, Column, Row As Integer
        'MSChart1.chartType = MSChart20Lib.VtChChartType.VtChChartType2dLine


        'MSChart20Lib.VtChChartType.VtChChartType2dLine

        'VtChChartType2dCombination()
        'MSChart20Lib.VtChChartType.VtChChartType2dXY()


        For i = 1 To 30
            arrData(i, 0) = i * i
        Next


        '        MSChart1.ChartData = arrData

    End Sub


    Private Sub Button1_Click_1(ByVal sender As System.Object, ByVal e As System.EventArgs)

    End Sub

    'Private Sub MSChart1_ChartSelected(ByVal sender As System.Object, ByVal e As AxMSChart20Lib._DMSChartEvents_ChartSelectedEvent)

    'End Sub

    'Private Sub MSChart1_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
    '   Call Dplot_Plot(s_Plot_File_Name(1), s_Plot_File_Type(1))

    'End Sub


    Private Sub Button1_Click_2(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(rootdata_directory & "test.dat", 0)
    End Sub
    Sub Dplot_Plot(ByRef Plot_File_name As String, ByRef Plot_Delta As Short)
        Dim d As DPLOT
        Dim i, ii, iii As Integer
        'Dim N As Integer
        Dim x() As Double
        Dim y() As Double
        Dim cmds As String
        Dim ret As Integer

        Dim Number_Measuring_Lights As Integer
        Dim total_points As Integer
        Dim j As String
        Dim Plot_Data_X() As Double
        Dim Plot_Data_Y() As Double
        Dim temp As Single = 0

        ' figure out how many columns


        FileOpen(1, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        j = LineInput(1)
        FileClose(1)
        'Threading.Thread.sleep(10)


        Number_Measuring_Lights = 1
        For i = 1 To Len(j)
            If Asc(Mid(j, i, 1)) = 9 Then  'counting the number of tabs
                Number_Measuring_Lights = Number_Measuring_Lights + 1
            End If
        Next i


        Number_Measuring_Lights = Number_Measuring_Lights / 4


        FileOpen(1, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)


        i = 0

        While EOF(1) = False
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            Input(1, temp)
            i = i + 1
        End While  'coiunt the total number of points, each has 4 data

        FileClose(1)
        'Threading.Thread.sleep(10)

        total_points = i

        ReDim Plot_Data_X(total_points)
        ReDim Plot_Data_Y(total_points)

        i = 0

        FileOpen(1, Plot_File_name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
        If Plot_Delta = 1 Then

            While EOF(1) = False
                Input(1, Plot_Data_X(i))
                Input(1, temp)
                Input(1, temp)
                Input(1, Plot_Data_Y(i))
                i = i + 1
            End While

        Else

            While EOF(1) = False
                Input(1, Plot_Data_X(i))
                Input(1, Plot_Data_Y(i))
                Input(1, temp)
                Input(1, temp)
                i = i + 1
            End While

        End If
        FileClose(1)
        'Threading.Thread.sleep(10)

        ' Plot sin(PI*x) and cos(PI*x) from x = 0 to 4.
        ' In this case we've used DATA_XYYY (one X array for one or more Y arrays).
        '
        ' Since the X values are evenly spaced we could also use
        ' DATA_DXY  - X has only 2 elements, DX and X0
        '
        ' And you can ALWAYS use
        ' DATA_XYXY - Each curve uses its own X array. This is the only
        '             option available if the curves have different X values or
        '             a different number of points.

        'N = total_points
        ReDim x(total_points)
        ReDim y(total_points)

        iii = 0
        For ii = 0 To Number_Measuring_Lights - 1

            For i = ii To (total_points - Number_Measuring_Lights) + ii Step Number_Measuring_Lights
                x(iii) = Plot_Data_X(i)
                y(iii) = Plot_Data_Y(i)
                iii = iii + 1
            Next i
        Next ii
        'For i = 0 To N - 1
        ' x(i) = Plot_Data_X(i) '(4.0# * i) / (N - 1)
        ' y(i) = Plot_Data_Y(i) 'System.Math.Sin(Math.PI * x(i))
        ' 'y(i + N) = Plot_Data_Y(i) * 2
        ' Next i

        d = New DPLOT
        d.Initialize()
        d.Version = DPLOT_DDE_VERSION
        d.DataFormat = DATA_XYXY
        d.MaxCurves = 16 ' Must be >= number of curves we plot
        d.MaxPoints = total_points * 2 ' Anything >= N will do
        d.NumCurves = Number_Measuring_Lights
        d.ScaleCode = SCALE_LINEARX_LINEARY
        d.LegendX = 0.05
        d.LegendY = 0.05
        For i = 0 To Number_Measuring_Lights - 1
            d.NP(i) = total_points / Number_Measuring_Lights
        Next i

        For i = 0 To Number_Measuring_Lights - 1
            d.LineType(i) = LINESTYLE_SOLID
            'd.NP(i) = total_points / Number_Measuring_Lights
        Next i

        '        d.LineType(0) = LINESTYLE_SOLID
        '        d.LineType(1) = LINESTYLE_LONGDASH
        d.Title1 = "Data sent to DPlot via DPLOTLIB.DLL"
        d.XAxis = "time (s)"
        d.YAxis = "y"

        ' The command string can be as complex as you like, limited to 32K characters.
        ' As an alternative, if your program will be producing many similar plots then
        ' you might prefer to create a preferences file, edit the file with a text
        ' editor to remove unwanted entries, and use the [GetPreferences("filename")]
        ' command.
        cmds = "" '"[ManualScale(0,-1.25,4,1.25)][TickInterval(1,0.5,0.25)]"
        cmds = cmds & "[Caption(" & Chr(34) & "DPLOTLIB XY Test" & Chr(34) & ")]"
        'cmds = cmds & "[Legend(1," & Chr(34) & "sin({\sp}x)" & Chr(34) & ")]"
        'cmds = cmds & "[Legend(2," & Chr(34) & "cos({\sp}x)" & Chr(34) & ")]"
        cmds = cmds & "[DocMaximize()][ClearEditFlag()]"
        'cmds = cmds & "Color(index,r,g,b)"

        '        ret = DPlot_Plot8(d, Plot_Data_X(0), Plot_Data_Y(0), cmds)
        ret = DPlot_Plot8(d, x(0), y(0), cmds)
        'ret = DPlot_AddData(d, Plot_Data_X(0), Plot_Data_Y(0), cmds)

        'List1.Items.Add(ret & " " & Number_Measuring_Lights & " " & total_points)
    End Sub

    Private Sub Button2_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)

    End Sub

    Private Sub MSChart2_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(s_Plot_File_Name(2), s_Plot_File_Type(2))
    End Sub

    Private Sub MSChart3_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(s_Plot_File_Name(3), s_Plot_File_Type(3))
    End Sub

    Private Sub MSChart4_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(s_Plot_File_Name(4), s_Plot_File_Type(4))
    End Sub

    Private Sub MSChart5_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(s_Plot_File_Name(5), s_Plot_File_Type(5))
    End Sub

    Private Sub MSChart6_DblClick(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call Dplot_Plot(s_Plot_File_Name(6), s_Plot_File_Type(6))
    End Sub

    Private Sub Label5_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Label5.Click

    End Sub

    Private Sub CreateGraph(ByVal zgc As ZedGraphControl, ByVal number_measuring_lights As Integer, ByVal total_points As Integer, ByVal plot_data_x As Object, ByVal plot_data_y As Object, ByVal plot_delta As Short, ByVal Graph_Name_l As Object, ByVal add_to As Integer, ByVal plot_window As Integer)
        Dim myPane As GraphPane = zgc.GraphPane
        Dim i, ii As Integer
        Dim x As Double
        Dim y As Double

        ' Set the titles and axis labels
        myPane.Title.Text = Graph_Name_l
        'Me.List1.Items.Add(Graph_Name_l)
        myPane.XAxis.Title.Text = "data point" '"time (s)"
        If plot_delta = 1 Then
            myPane.YAxis.Title.Text = ".delta.A"
        Else
            myPane.YAxis.Title.Text = "V"
        End If

        If add_to = False Then
            myPane.CurveList.Clear()
            ccc(plot_window) = 0
        End If


        ' Make up some data points from the Sine function


        ', y As Double

        For ii = 1 To number_measuring_lights
            ccc(plot_window) = ccc(plot_window) + 1
            Dim list1 As New PointPairList()

            For i = 0 To total_points - 1
                If plot_data_point = True Then
                    x = i ' plot_data_x(ii, i)
                Else
                    x = plot_data_x(ii, i)
                End If


                y = plot_data_y(ii, i)
                list1.Add(x, y)
            Next i

            Select Case (ccc(plot_window) And 7)
                Case 1
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Red, SymbolType.None)
                Case 2
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Orange, SymbolType.None)
                Case 3
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Yellow, SymbolType.None)
                Case 4
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Green, SymbolType.None)
                Case 5
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Blue, SymbolType.None)
                Case 6
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.Violet, SymbolType.None)
                Case 7
                    Dim myCurve1 As LineItem = myPane.AddCurve("T:" & Str(Current_Trace), list1, Color.White, SymbolType.None)

                    'myCurve1.Symbol.Fill = New Fill(Color.White)
            End Select
        Next
        ' Generate a blue curve with circle symbols, and "My Curve 2" in the legend

        ' Fill the axis background with a color gradient
        myPane.Chart.Fill = New Fill(Color.Black, Color.Green, 45.0F)

        ' Fill the pane background with a color gradient
        myPane.Fill = New Fill(Color.White, Color.FromArgb(255, 255, 155), 45.0F)

        ' Calculate the Axis Scale Ranges
        zgc.AxisChange()

        zgc.Refresh()

    End Sub



    Private Sub zg1_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles zg1.Load

    End Sub



    Private Sub Button7_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button7.Click
        'Exit Sub
        'Dim ulstat As Integer
        'Dim direction As Integer = 0
        'Dim enable As Integer = 0
        'Dim p_v As Short
        'Dim p_o, P_O_A As Single
        'Dim target As Single
        'Dim i, ii As Integer
        'target = Val(TextBox1.Text)
        'p_o = 1
        'P_O_A = 1

        'ulstat = cbAIn(BoardNum, 7, BIP5VOLTS, p_v)
        'ulstat = cbToEngUnits(BoardNum, BIP5VOLTS, p_v, p_o)
        'Label1.Text = p_o

        'Exit Sub
        'Halt_Script = False


        'If target > Current_Wheel_Position Then
        '    direction = 1
        'Else
        '    direction = 0
        'End If

        'For i = 1 To Math.Abs(Current_Wheel_Position - target)
        '    enable = 1
        '    D_Out_Mask = 1 * enable + 2 * direction
        '    ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
        '    enable = 0
        '    D_Out_Mask = 1 * enable + 2 * direction

        '    ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)

        '    ulstat = cbAIn(BoardNum, 7, BIP5VOLTS, p_v)
        '    ulstat = cbToEngUnits(BoardNum, BIP5VOLTS, p_v, p_o)

        '    'Application.DoEvents()
        '    For ii = 1 To 100000
        '    Next ii

        'Next i
        'Current_Wheel_Position = target
        'Label1.Text = Current_Wheel_Position


        'Exit Sub

        'While P_O_A > target + 0.001 Or P_O_A < target - 0.001
        '    If Halt_Script = True Then Exit While
        '    If p_o > target Then
        '        direction = 1
        '    Else
        '        direction = 0
        '    End If

        '    enable = 1
        '    D_Out_Mask = 1 * enable + 2 * direction
        '    ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
        '    enable = 0
        '    D_Out_Mask = 1 * enable + 2 * direction
        '    ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
        '    'ulstat=
        '    P_O_A = 0
        '    For i = 1 To 8
        '        ulstat = cbAIn(BoardNum, 7, BIP5VOLTS, p_v)
        '        ulstat = cbToEngUnits(BoardNum, BIP5VOLTS, p_v, p_o)
        '        P_O_A = P_O_A + p_o
        '    Next i
        '    P_O_A = P_O_A / 8
        '    Label1.Text = P_O_A
        '    Application.DoEvents()

        'End While

    End Sub

    Private Sub TextBox1_TextChanged(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TextBox1.TextChanged

    End Sub

    Private Sub Label1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Label1.Click

    End Sub
    Private Sub move_wheel(ByVal set_position As Single, ByVal av_it As Integer)
        Dim ulstat, direction, enable As Integer
        Dim p_v As Short
        Dim i As Integer
        Dim p_o, P_O_A As Single
        Dim target As Single
        target = set_position 'Val(TextBox1.Text)
        ulstat = cbAIn(BoardNum, 7, BIP10VOLTS, p_v)
        ulstat = cbToEngUnits(BoardNum, BIP10VOLTS, p_v, p_o)
        Label1.Text = p_o
        Halt_Script = False
        While P_O_A > target + 0.001 Or P_O_A < target - 0.001
            If Halt_Script = True Then Exit While
            If p_o > target Then
                direction = 1
            Else
                direction = 0
            End If

            enable = 1
            D_Out_Mask = 1 * enable + 2 * direction
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
            enable = 0
            D_Out_Mask = 1 * enable + 2 * direction
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
            'ulstat=
            P_O_A = 0
            For i = 1 To av_it
                ulstat = cbAIn(BoardNum, 7, BIP10VOLTS, p_v)
                ulstat = cbToEngUnits(BoardNum, BIP10VOLTS, p_v, p_o)
                P_O_A = P_O_A + p_o
            Next i
            P_O_A = P_O_A / av_it

            Application.DoEvents()

        End While
        Label1.Text = P_O_A


    End Sub

    Private Sub HomewheelToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Call move_wheel(1.51, 1)
        Call move_wheel(1.51, 32)
        Current_Wheel_Position = 0
    End Sub
    Private Sub set_wheel(ByVal target As Integer)
        'On Error Resume Next
        Dim receive_data As String = ""
        Label1.Text = "moving"
        System.Windows.Forms.Application.DoEvents()
        Using com1 As IO.Ports.SerialPort = My.Computer.Ports.OpenSerialPort("COM1")

            com1.BaudRate = 19200
            com1.DataBits = 8
            com1.StopBits = 1
            com1.Parity = IO.Ports.Parity.None
            com1.ReadTimeout = 20000

            com1.WriteLine("WGOTO" & target)


            'While receive_data = "" And Halt_Script = False
            receive_data = com1.ReadLine()
            System.Windows.Forms.Application.DoEvents()
            'End While
            Label1.Text = receive_data
            If receive_data = "*" Then
                Label1.Text = "OK"
            Else
                Label1.Text = "TIMED OUT"
            End If

        End Using

    End Sub
    Private Sub home_wheel()
        'On Error Resume Next
        Using com1 As IO.Ports.SerialPort = My.Computer.Ports.OpenSerialPort("COM1")
            Dim receive_data As String = ""
            com1.BaudRate = 19200
            com1.DataBits = 8
            com1.StopBits = 1
            com1.Parity = IO.Ports.Parity.None
            com1.ReadTimeout = 20000

            com1.WriteLine("WSMODE")


            'While receive_data = "" And Halt_Script = False
            receive_data = com1.ReadLine()
            System.Windows.Forms.Application.DoEvents()
            'End While
            Label1.Text = receive_data

            If receive_data <> "!" Then
                Label1.Text = "TIMED OUT"
                Exit Sub
            End If


            com1.WriteLine("WHOME")


            'While receive_data = "" And Halt_Script = False
            receive_data = com1.ReadLine()
            System.Windows.Forms.Application.DoEvents()
            'End While

            Label1.Text = receive_data


            com1.Close()
        End Using
    End Sub
    Private Sub Simmer(ByVal number_pulses As Integer, ByVal q_switch As String, ByVal simmer_q_switch As Short)
        'Dim ulstat As Short

        'ulstat = cbDOut(BoardNum, AUXPORT, 0)

        '        Trace_Running = True

        'Dim Time_In_S As Single
        Dim PP_Line As Short = 0
        '       Dim i, j As Object
        Dim ii As Integer = 0
        Dim k As Short = 0
        Dim AD_Trigger As Short = 1

        'Dim Loop_Return_Index As Short

        Dim ss As Short = 0





        points = 0


        ' simmer protocol
        cblaster.CreateNewProtocol(0, 1)
        cblaster.AddNextBitmask(0, 0, cblaster.makeBitmask(0, 0, 0, 0, 0, 0, 0, 0, 1, 0), Protocol.TimeInSeconds("5u")) ' bit22 on
        cblaster.AddNextBitmask(0, 0, cblaster.makeBitmask(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), Protocol.TimeInSeconds(q_switch)) ' bit22 off
        cblaster.AddNextBitmask(0, 0, cblaster.makeBitmask(0, 0, 0, 0, 0, 1, 0, 0, 1, 0), Protocol.TimeInSeconds("10u")) 'bit22 on again, now with farred bit (bit19)
        cblaster.AddNextBitmask(0, 0, cblaster.makeBitmask(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), Protocol.TimeInSeconds("100m")) ' off again
        cblaster.SetPulseNumber(0, 0, number_pulses)
        cblaster.SetPostProtocolBitmask(0, 0)
        If Not cblaster.ExecuteProtocol(0, (Protocol.TimeInSeconds("5u") + Protocol.TimeInSeconds(q_switch) + _
                                Protocol.TimeInSeconds("10u") + Protocol.TimeInSeconds("100m")) * number_pulses * 2) Then
            List1.Items.Add("SIMMER failed")
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1
            System.Windows.Forms.MessageBox.Show("Error occured while programming the ChrisBlaster", "FPGA Communication Error")
            System.Windows.Forms.Application.DoEvents()
        End If




        List1.Items.Add(" SIMMER")
        'Scrolling text box fix: make the box jump to the bottom after adding an item
        List1.TopIndex = List1.Items.Count - 1
        System.Windows.Forms.Application.DoEvents()



    End Sub

    Private Sub Button8_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button8.Click
        Call Simmer(1000, (Str(TrackBar1.Value) & "u"), 0)
    End Sub

    Private Sub TrackBar1_Scroll(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TrackBar1.Scroll
        Label6.Text = TrackBar1.Value
    End Sub

    Private Sub Button9_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button9.Click
        Call Simmer(1000, (Str(TrackBar1.Value) & "u"), 1)
    End Sub

    Private Sub Button10_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button10.Click
        Dim test As String
        test = Trim("IDEA" & DateString & TimeOfDay)


        List1.Items.Add(test)

        Exit Sub


        'Dim x As String
        'Dim i As Integer = 0
        'Dim mscomm1 as com


        'x = Chr(255) & Chr(0) & Chr(100)

        'Using com1 As IO.Ports.SerialPort = My.Computer.Ports.OpenSerialPort("COM1")
        Using com1 As IO.Ports.SerialPort = My.Computer.Ports.OpenSerialPort("COM1")
            'My.Computer.Ports.OpenSerialPort("COM1")
            'com1.Encoding
            'com1.DtrEnable = True
            com1.Handshake = IO.Ports.Handshake.None
            com1.NewLine = 0

            com1.DtrEnable = False
            com1.Parity = IO.Ports.Parity.None
            'com1.

            com1.BaudRate = 2400
            com1.DataBits = 8
            com1.StopBits = IO.Ports.StopBits.One
            'com1.Encoding.ASCII(
            com1.WriteTimeout = 100000


            'com1.Write(x)
            'com1.Write (x)
            'com1.Write(Chr(255) & Chr(0) & Chr(22))
            'com1.
            'For i = 0 To 255
            'com1.Write(Chr(255) & Chr(0) & Chr(i))
            'Next i
            'Label1.Text = Asc((Chr(255))) & Asc(Chr(0)) & Asc(Chr(10))

            com1.Write(Chr(255))
            Call Hold_on_There(0.2)
            com1.Write(Chr(0))
            Call Hold_on_There(0.2)
            com1.Write(Chr(0))

            '           com1.WriteLine(Chr(255))
            '           com1.WriteLine(Chr(0))
            '          com1.WriteLine(Chr(25))

            'Call Hold_on_There(0.2)
            'Next

            '            Label1.Text = Asc((Chr(255))) & Asc(Chr(0)) & Asc(Chr(100))


            'com1.WriteLine(Chr(255))
            'com1.Write(Chr(255))
            'com1.Write(Chr(255))
            'com1.Write(Chr(255))
            'com1.Write(Chr(255))

            'com1.WriteLine(Chr(1))

            'com1.WriteLine(Chr(255))


        End Using

        Exit Sub

    End Sub

    Private Sub OnToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)

    End Sub

    Private Sub LaserToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        If view_laser = True Then
            Button8.Visible = False
            Button9.Visible = False
            TrackBar1.Visible = False
            Label6.Visible = False

            view_laser = False
        Else
            view_laser = True
            Button8.Visible = True
            Button9.Visible = True
            TrackBar1.Visible = True
            Label6.Visible = True
        End If

    End Sub

    Private Sub Button11_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        'Call set_servo(1, Val(TextBox2.Text))
    End Sub
    Private Sub set_xe_intensity(ByVal inv As Single)
        If inv > 10 Then inv = 9.9
        If inv < 0 Then inv = 0


        Dim ttt As Long
        Dim outv, t As Short
        ttt = inv * (65535) / 10
        If ttt < 32768 Then
            outv = ttt
        Else
            outv = (ttt - 65535)

        End If
        'Label1.Text = "*" & outv
        '~~~

        t = cbAOut(BoardNum, 0, 1, outv)
        System.Windows.Forms.Application.DoEvents()
    End Sub

    Sub set_servo(ByVal servo_num As Short, ByVal inv As Single)
        '        Label7.Text = (servo_num & " " & inv)

        Exit Sub

        'Dim inv As Single
        Dim ttt As Long
        Dim outv As Short = 0
        Dim t As Short = 0
        'inv = Val(TextBox1.Text)
        ttt = (1000 * inv / 180) + 250
        'inv * (65535) / 10
        'If ttt < 32768 Then
        ' outv = ttt
        ' Else
        ' outv = (ttt - 65535)

        'End If
        'Label2.Text = outv

        't = cbAOut(1, 0, 1, outv)


        t = cbC9513Init(1, 1, 2, FREQ2, CBDISABLED, CBDISABLED, CBDISABLED)

        t = cbC9513Config(1, servo_num, NOGATE, POSITIVEEDGE, FREQ2, CBDISABLED, LOADANDHOLDREG, RECYCLE, CBDISABLED, COUNTDOWN, TOGGLEONTC)
        t = cbCLoad(1, LOADREG1, 12500)
        't = cbCLoad(2, HOLDREG1, (1000 - (1000 * inv)))
        t = cbCLoad(1, HOLDREG1, (ttt))
        't = cbtimerout
    End Sub

    Private Sub FilterwheelToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs)
        Form3.Show()

    End Sub

    Private Sub TakeFluorBaselineToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TakeFluorBaselineToolStripMenuItem.Click
        Dim ulstat As Integer
        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)
        If fluorescence_shutter_status = 0 Then
            Call fluorescence_shutter("1")

        Else
            Call fluorescence_shutter("0")

        End If

    End Sub

    Private Sub fluorescence_shutter(ByVal set_at As String)
        Dim ulstat As Integer
        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)
        If set_at = "open" Or set_at = "1" Then
            fluorescence_shutter_status = 1

            D_Out_Mask = (D_Out_Mask Or 2)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)

            Call Hold_on_There(Shutter_delay)
            D_Out_Mask = (D_Out_Mask And &HFD)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
            List1.Items.Add(fluorescence_shutter_status)
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1

        ElseIf set_at = "close" Or set_at = "0" Then
            fluorescence_shutter_status = 0
            D_Out_Mask = (D_Out_Mask Or 1)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)

            Call Hold_on_There(Shutter_delay)
            D_Out_Mask = (D_Out_Mask And &HFE)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)

            List1.Items.Add(fluorescence_shutter_status)
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1

        End If
    End Sub

    Private Sub stir_control(ByVal on_off As Short)
        Dim ulstat As Integer
        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)
        If on_off > 0 Then
            stir_status = 1

            D_Out_Mask = (D_Out_Mask Or 4)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)



        ElseIf on_off = 0 Then
            stir_status = 0
            D_Out_Mask = (D_Out_Mask And &HB)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)


        End If
    End Sub

    Public Sub TakenoteToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TakenoteToolStripMenuItem.Click
        Note_File_Name = Trim(Base_File_Name & "_notes.txt")
        take_a_note.Show()

    End Sub

    Private Sub NotesToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles NotesToolStripMenuItem.Click

    End Sub

    Private Sub SetbasefileToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles SetbasefileToolStripMenuItem.Click
        Call set_base_file_name()
    End Sub

    Private Sub RecordeventsToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles RecordeventsToolStripMenuItem.Click

    End Sub

    Private Sub ONToolStripMenuItem_Click_1(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles ONToolStripMenuItem.Click
        take_notes = True

    End Sub

    Private Sub OFFToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles OFFToolStripMenuItem.Click
        take_notes = False
    End Sub
    'Public Sub do_auto_gain(ByRef current_protocol As Integer, ByRef M_Reference_Gain(,) As Short, ByRef L_Measuring_Interval(,,) As String, _
    '                        ByRef m_measuring_light(,) As Short, ByRef m_choke_actinic(,) As Short, ByRef m_number_measuring_lights() As Short, _
    '                        ByRef m_number_loops() As Short, ByRef M_Intensity(,) As Short, ByRef M_Number_Pulses(,) As Short, ByRef M_Take_Data(,) _
    '                        As Short, ByRef Measuring_Pulse_Duration(,) As String, ByRef Xe_Flash(,) As Short, ByRef Gain() As Integer, _
    '                        ByRef In_Channel() As Short, ByRef Ref_Channel_Number As Short, ByRef M_Far_Red(,) As Short, ByRef M_Blue_Actinic(,) As Short, ByRef gain_slop As Single, _
    '                        ByRef q_switch(,) As Short, ByRef Pre_Pulse(,) As Short, ByRef Pre_Pulse_Time(,) As String, ByRef Pre_Delay(,) As String)

    '    'do_auto_gain(M_Number_Measuring_Lights(Current_Protocol), 
    '    Dim ii, iiii As Integer
    '    Dim Max_Volts As Single
    '    Dim gain_temp, gain_set_temp As Single
    '    Dim gain_data(10) As Single
    '    'gain_data_list.Items.Clear()
    '    For iiii = 1 To m_number_measuring_lights(current_protocol)
    '        m_number_loops(0) = 1
    '        For ii = 1 To m_number_loops(0)
    '            M_Number_Pulses(0, ii) = 1
    '            M_Intensity(0, ii) = M_Intensity(current_protocol, 1)
    '        Next ii
    '        M_Take_Data(0, 1) = 1
    '        'm_number_measuring_lights(0) = 8

    '        For ii = 1 To m_number_measuring_lights(0)
    '            m_measuring_light(0, ii) = m_measuring_light(current_protocol, iiii)  'set the measuring light to each of those used in the protocol
    '        Next ii
    '        For ii = 1 To m_number_measuring_lights(0)
    '            m_choke_actinic(0, ii) = m_choke_actinic(current_protocol, iiii)  'set the measuring light to each of those used in the protocol
    '        Next ii

    '        '           For ii = 1 To m_number_measuring_lights(0)
    '        ' M_Measuring_Interval(0, ii) = "10m"
    '        'Next ii
    '        For ii = 1 To m_number_measuring_lights(0)
    '            L_Measuring_Interval(0, ii, 1) = "10m"
    '        Next ii
    '        'For ii = 1 To m_number_measuring_lights(0)
    '        '    M_Primary_Gain(0, ii) = ii - 1
    '        'Next ii
    '        'For ii = 1 To m_number_measuring_lights(0)
    '        '    M_Reference_Gain(0, ii) = ii - 1
    '        'Next ii
    '        For ii = 1 To m_number_measuring_lights(0)
    '            Measuring_Pulse_Duration(0, ii) = "10u" 'Str$(ii * 20) & "u"
    '        Next ii

    '        Gain(0) = Gain(current_protocol)
    '        In_Channel(0) = In_Channel(current_protocol)
    '        'Xe_Flash, M_Far_Red, M_Blue_Actinic
    '        For ii = 1 To m_number_loops(0)
    '            Xe_Flash(0, ii) = Xe_Flash(current_protocol, 1)
    '            M_Far_Red(0, ii) = M_Far_Red(current_protocol, 1)
    '            M_Blue_Actinic(0, ii) = M_Blue_Actinic(current_protocol, 1)
    '        Next ii
    '        '                    Call Multi_Trace(0, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, M_Measuring_Light, M_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, Gain(Current_Protocol), In_Channel(Current_Protocol), points, data_time)
    '        Call Multi_Trace(0, m_number_loops, M_Intensity, M_Number_Pulses, m_number_measuring_lights, m_measuring_light, m_choke_actinic, _
    '                         L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, Measuring_Pulse_Duration, _
    '                         Gain, In_Channel, Ref_Channel_Number, _
    '                         points, data_time, Xe_Flash, q_switch, M_Far_Red, _
    '                         M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, M_Take_Data)
    '        'Call Dirk_Save(data_time, points, "c:\gain.dat", File_Replace, 1, 1, 0, m_number_measuring_lights(0), 1, 1, In_Channel(0))

    '        'Call Plot_File("c:\gain.dat", 1, (0), 0, 1)
    '        'Call Plot_File("c:\gain.dat", 2, (2), 0, 1)

    '        'For ii = 0 To m_number_measuring_lights(0) - 1
    '        'gain_data(ii) = Data_Volts(ii, In_Channel(current_protocol))
    '        'gain_data_list.Items.Add("gain data: " & M_Primary_Gain(0, ii + 1) & ": " & gain_data(ii))
    '        'Next ii

    '        'For ii = 0 To m_number_measuring_lights(0) - 1
    '        'gain_data(ii) = Data_Volts(ii, 2) 'reference channel
    '        'gain_data_list.Items.Add("gain ref: " & M_Reference_Gain(0, ii + 1) & ": " & gain_data(ii))
    '        ' Next ii
    '        'BIP10VOLTS, bip5volts, bip2pt5volts, bip1volts

    '        If Gain(0) = BIP10VOLTS Then
    '            Max_Volts = 10 * gain_slop
    '        ElseIf Gain(0) = BIP5VOLTS Then
    '            Max_Volts = 5 * gain_slop
    '        ElseIf Gain(0) = BIP2VOLTS Then
    '            Max_Volts = 2 * gain_slop
    '        ElseIf Gain(0) = BIP1VOLTS Then
    '            Max_Volts = 1 * gain_slop
    '        End If
    '        gain_temp = 0

    '        For ii = 0 To m_number_measuring_lights(0) - 1  'determine the best gain for the data channel
    '            gain_data(ii) = Data_Volts(ii, In_Channel(current_protocol))
    '            If Math.Abs(gain_data(ii)) > Math.Abs(gain_temp) And Math.Abs(gain_data(ii)) < Max_Volts Then
    '                'If Math.Abs(gain_data(ii)) < Max_Volts Then

    '                gain_temp = gain_data(ii)
    '                gain_set_temp = ii
    '            End If
    '        Next ii

    '        'gain_data_list.Items.Add("best data: " & gain_temp & ": " & gain_set_temp)
    '        M_Primary_Gain(current_protocol, iiii) = gain_set_temp
    '        ProgressBar1.Value = M_Primary_Gain(current_protocol, iiii)
    '        sample_gain_label.Text = M_Primary_Gain(current_protocol, iiii)
    '        gain_temp = 0
    '        For ii = 0 To m_number_measuring_lights(0) - 1  'determine the best gain for the ref channel
    '            gain_data(ii) = Data_Volts(ii, 2)
    '            If Math.Abs(gain_data(ii)) > Math.Abs(gain_temp) And Math.Abs(gain_data(ii)) < Max_Volts Then
    '                'If Math.Abs(gain_data(ii)) < Max_Volts Then
    '                gain_temp = gain_data(ii)
    '                gain_set_temp = ii
    '            End If
    '        Next ii

    '        'gain_data_list.Items.Add("best ref: " & gain_temp & ": " & gain_set_temp)
    '        M_Reference_Gain(current_protocol, iiii) = gain_set_temp
    '        ProgressBar2.Value = M_Reference_Gain(current_protocol, iiii)
    '        reference_gain_label.Text = M_Reference_Gain(current_protocol, iiii)

    '    Next iiii

    'End Sub
    'Public Sub do_auto_gain_ref(ByVal light As Integer)
    '    '        Call do_auto_gain(Current_Protocol, M_Reference_Gain, L_Measuring_Interval, M_Measuring_Interval, M_Measuring_Light, M_Number_Measuring_Lights, M_Number_Loops, M_Intensity, M_Number_Pulses, M_Take_Data, Measuring_Pulse_Duration, Xe_Flash, Gain, In_Channel, M_Far_Red, M_Blue_Actinic, gain_slop, q_switch, Pre_Pulse, Pre_Pulse_Time, Pre_Delay)
    '    script_file = rootdata_directory & "gain_test.txt"
    '    'M_Measuring_Light(1, 1) = light
    '    Call run_the_script()
    '    script_file = ""
    'End Sub

    Private Sub OnToolStripMenuItem1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles OnToolStripMenuItem1.Click
        Call stir_control(1)
    End Sub

    Private Sub OffToolStripMenuItem1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles OffToolStripMenuItem1.Click
        Call stir_control(0)
    End Sub

    Private Sub StirerToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles StirerToolStripMenuItem.Click

    End Sub
    ' This method is designed to produce more infomrative script-reading error messages
    Protected Function ShowScriptErrorMessage(ByVal linenumber As Integer, ByRef lines() As String, ByRef message As String) As Boolean
        Dim margin As Integer = 5
        Dim start As Integer = linenumber - margin
        Dim ans As DialogResult
        If start < 0 Then
            start = 0
        End If
        Dim code As String = vbNewLine
        While (start < (linenumber + margin) And start < lines.Length)
            code = code & lines(start) & vbNewLine
            start += 1
        End While

        ans = System.Windows.Forms.MessageBox.Show("Error reading script file:" & vbNewLine & "Command: """ & lines(linenumber) & """ (read #" & linenumber & ")" & vbNewLine & _
                              message & vbNewLine & "Error found in code block" & vbNewLine & code & vbNewLine & _
                              "Continue?", "READ ERROR", MessageBoxButtons.YesNo)
        If ans = Windows.Forms.DialogResult.No Then
            Return False
        Else
            Return True
        End If
    End Function
    Protected Function ShowScriptFatalErrorMessage(ByVal linenumber As Integer, ByRef lines() As String, ByRef message As String) As Boolean
        Dim margin As Integer = 5
        Dim start As Integer = linenumber - margin
        Dim ans As DialogResult
        If start < 0 Then
            start = 0
        End If
        Dim code As String = vbNewLine
        While (start < (linenumber + margin) And start < lines.Length)
            code = code & lines(start) & vbNewLine
            start += 1
        End While

        ans = System.Windows.Forms.MessageBox.Show("Error reading script file:" & vbNewLine & "Command: """ & lines(linenumber) & """ (read #" & linenumber & ")" & vbNewLine & _
                              message & vbNewLine & "Error found in code block" & vbNewLine & code & vbNewLine, _
                                "READ ERROR", MessageBoxButtons.OK)
        If ans = Windows.Forms.DialogResult.No Then
            Return False
        Else
            Return True
        End If
    End Function

    'This is the code generated for the Test_Script menu item
    Private Sub TestScriptToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles TestScriptToolStripMenuItem.Click
        List1.Items.Add("TESTING SCRIPT FILE")
        List1.TopIndex = List1.Items.Count - 1
        Me.BackColor = Color.Yellow
        Call run_the_script(False)
        Me.BackColor = Color.Gray
        List1.Items.Add("TEST COMPLETE")
        List1.TopIndex = List1.Items.Count - 1
    End Sub

    Private Sub Button11_Click_1(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button11.Click
        List1.Items.Add("TESTING SCRIPT FILE")
        List1.TopIndex = List1.Items.Count - 1
        Me.BackColor = Color.Yellow
        Call run_the_script(False)
        Me.BackColor = Color.Gray
        List1.Items.Add("TEST COMPLETE")
        List1.TopIndex = List1.Items.Count - 1
    End Sub

    ' subroutine for testing the saturation light (0.5 second flash)
    Public Sub testSaturationLight()
        ' check if the test protocol has been created
        If IsNothing(saturation_test) Then
            ' create the test protocol
            saturation_test = New Protocol()
            saturation_test.Number_Loops = 2
            saturation_test.Number_Measuring_Lights = 1
            saturation_test.Measuring_Light(0) = 0
            For L As Integer = 0 To saturation_test.Number_Loops - 1
                saturation_test.Saturating_Pulse(L) = 1
                saturation_test.Take_Data(L) = 0
                saturation_test.Measuring_Interval(0, L) = "0.01" ' units in 100ths of a second
                saturation_test.Number_Pulses(L) = 50
                saturation_test.Intensity(L) = 0
                saturation_test.Far_Red(L) = 0
                saturation_test.Blue_Actinic(L) = 0
            Next
        End If
        Dim ignore(-1) As Double
        cblaster.ProgramProtocol(0, saturation_test, ignore)
        cblaster.StartRunningProtocol(0)
    End Sub

    Private Sub SetAutoRunFileToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles SetAutoRunFileToolStripMenuItem.Click
        If auto_run = False Then
            auto_run = True
            Call set_auto_file_name()
        Else
            auto_run = False
        End If
    End Sub

    Private Sub Timer1_Tick(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Timer1.Tick
        Dim temp As String = ""
        Dim temp2 As String = ControlChars.CrLf
        Dim temp3 As String = ControlChars.CrLf

        On Error GoTo skip

        If auto_run = True Then
            FileOpen(6, Auto_File_Name, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
            temp = LineInput(6)

            FileClose(6)
            If temp <> "" Then
                Timer1.Enabled = False
                List1.Items.Add(Auto_File_Name)

                script_file = Auto_File_Name

                Script_Label.Text = script_file
                run_another_script = True
                'script_file = Script_Label.Text
                'Me.Text = "Kinetic Spectrophotometer: RUNNING"
                Me.BackColor = Color.Green
                While run_another_script = True
                    Call run_the_script()
                End While
                'Me.Text = "Kinetic Spectrophotometer: IDLE"
                Me.BackColor = Color.Gray

                FileOpen(1, script_file, OpenMode.Input, OpenAccess.Default, OpenShare.Shared)
                'temp3 = ControlChars.CrLf
                While EOF(1) = False
                    temp2 = LineInput(1)
                    temp3 = temp3 & temp2 & ControlChars.CrLf

                End While
                FileClose(1)

                FileOpen(7, script_file, OpenMode.Output, OpenAccess.Default, OpenShare.Shared)
                Write(7, temp3)
                FileClose(7)

            End If

            ' List1.Items.Add(temp)
            Exit Sub
        End If
skip:
        '        List1.Items.Add("skip " & Auto_File_Name)
    End Sub

    Private Sub Cd1ToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Cd1ToolStripMenuItem.Click
        'Imports Microsoft.VisualBasic.FileIO.TextFieldParser
        'Imports System.IO

        Dim special_function_file As String = ""
        CommonDialog1Open.Title = "enter file name (& location) for SPECIAL"
        CommonDialog1Save.Title = "enter file name (& location) for SPECIAL"
        CommonDialog1Open.FileName = "*.txt"
        CommonDialog1Save.FileName = "*.txt"
        CommonDialog1Open.Title = "SPECIAL DATA ANALYSIS FUNCTION"
        CommonDialog1Save.Title = "SPECIAL DATA ANALYSIS FUNCTION"
        CommonDialog1Save.ShowDialog()
        CommonDialog1Open.FileName = CommonDialog1Save.FileName


        If CommonDialog1Open.FileName = "" Then
            'Auto_File_Name = "temp"

        Else
            special_function_file = CommonDialog1Open.FileName
            'Base_File_Name = CommonDialog1Open.FileName
            List1.Items.Add("& " & special_function_file)
            'Scrolling text box fix: make the box jump to the bottom after adding an item
            List1.TopIndex = List1.Items.Count - 1
            If File.Exists(special_function_file) Then

                Dim afile As FileIO.TextFieldParser = New FileIO.TextFieldParser(special_function_file)
                Dim CurrentRecord(10) As String ' this array will hold each line of data
                'Tell the parser we are using delimited text file
                afile.TextFieldType = FileIO.FieldType.Delimited
                'Tell Parser how we will delimit
                afile.Delimiters = New String() {vbTab}
                'For this example field will be enclosed in quotes i.e. "MyEntry","MyEntry2"
                afile.HasFieldsEnclosedInQuotes = False

                Do While Not afile.EndOfData
                    Try
                        CurrentRecord = afile.ReadFields
                        '.ReadFields()
                        'Now we have our line in an array CurrentRecord 
                        'Do what we wish with it
                        List1.Items.Add(CurrentRecord)
                    Catch ex As FileIO.MalformedLineException
                        MsgBox("Malformed Text")
                        Exit Do
                    End Try
                Loop

            Else
                MsgBox("File does not exist")
            End If

        End If
        Exit Sub
err_it3:

        Err.Clear()
        script_file = ""
        Script_Label.Text = script_file
        'Resume get_out
get_out:
    End Sub

    Private Sub CustomdataToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles CustomdataToolStripMenuItem.Click

    End Sub

    Private Sub MatToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MatToolStripMenuItem.Click
        Dim receive_data As String

        'Bug_Out = ""
        Label1.Text = "moving"
        System.Windows.Forms.Application.DoEvents()
        Using com3 As IO.Ports.SerialPort = My.Computer.Ports.OpenSerialPort("COM3")

            com3.BaudRate = 9600
            com3.DataBits = 8
            com3.StopBits = 0
            com3.Parity = IO.Ports.Parity.None
            com3.ReadTimeout = 10000

            '            com1.WriteLine("WSMODE")


            com3.WriteLine(TextBox2.Text)

            System.Windows.Forms.Application.DoEvents()
            receive_data = com3.ReadLine()

            Label1.Text = receive_data
            If receive_data = "*" Then
                List1.Items.Add("OK")

            Else
                List1.Items.Add("TIMED OUT")
            End If
        End Using
    End Sub





    Private Sub Button12_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button12.Click
        Dim i, ulstat As Integer
        ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)

        For i = 1 To Val(TextBox3.Text)

            Hold_on_There(0.001)

            D_Out_Mask = (D_Out_Mask Or 64)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
            Label12.Text = i & "  " & D_Out_Mask

            Hold_on_There(0.001)

            D_Out_Mask = (D_Out_Mask And 191)
            ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)


        Next
    End Sub

    Private Sub Set_Lambda(ByVal target As Single)
        Exit Sub
        Dim ulstat As Integer
        Dim direction As Integer = 0
        Dim enable As Integer = 0
        Dim ans As String
        Dim steps_to_move As Integer
        Dim nanometers_actually_moved As Single
        If spec_mode = True Then
            ulstat = cbDConfigPort(mono_BoardNum, AUXPORT, DIGITALOUT)
        End If
        'Dim target As Single
        Dim i As Integer

        If target < 300 Or target > 750 Then
            ans = System.Windows.Forms.MessageBox.Show("Enter correct value for set lambda", "READ ERROR", MessageBoxButtons.OK)
            Exit Sub
        End If
        'current_lambda = Val(LambdaText.Text)
        If current_lambda < 300 Or current_lambda > 730 Then
            ans = System.Windows.Forms.MessageBox.Show("Enter correct value for monochromator lambda", "READ ERROR", MessageBoxButtons.OK)
            Exit Sub
        End If
        'ulstat = cbDConfigPort(mono_BoardNum, AUXPORT, DIGITALOUT)

        If target > current_lambda Then
            direction = 0
        Else
            direction = 1
        End If
        'Label4.Text = target
        System.Windows.Forms.Application.DoEvents()
        steps_to_move = (current_lambda - target) * steps_per_nanometer

        For i = 1 To Math.Abs(steps_to_move)
            enable = 1
            mono_D_Out_Mask = &H80 * enable + &H40 * direction
            ulstat = cbDOut(mono_BoardNum, AUXPORT, mono_D_Out_Mask)
            Hold_on_There(0.007)
            enable = 0
            mono_D_Out_Mask = &H80 * enable + &H40 * direction
            ulstat = cbDOut(mono_BoardNum, AUXPORT, mono_D_Out_Mask)

            Hold_on_There(0.007)

        Next i
        nanometers_actually_moved = steps_to_move / steps_per_nanometer
        Label4.Text = nanometers_actually_moved
        current_lambda = current_lambda - nanometers_actually_moved
        LambdaText.Text = current_lambda
        'My.Settings.lambda = current_lambda
        'My.Settings.Save()

    End Sub

    Private Sub Button13_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button13.Click
        'current_lambda = Val(LambdaText.Text)
        'Set_Lambda(Val(lambdasettext.Text))
    End Sub

    Private Sub GetSerialPortNames()
        'For Each sport As String In My.Computer.Ports.SerialPortNames
        ' cmbPort.Items.Add(sport)
        ' Next
    End Sub


    Private Sub Button14_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button14.Click
        'Set_Geoffs_detector(mdsteps.Text, Val(mddir.Text))
    End Sub

    ' Sub move_geoff_linear(ByVal move_cmd As String) ', ByVal holdoff As Integer)
    '     Dim direction As Integer
    '     Dim holdoff As Integer = 10000
    '     Dim number_of_steps_to_move As Integer = 0
    '     If move_cmd = "home" Then

    '     Else
    '         number_of_steps_to_move = Int(move_cmd)
    '     End If

    '     If number_of_steps_to_move > 1 Then
    '         direction = 1
    '     Else
    '         direction = 0
    '     End If


    '     Dim limit_switch As Single = 5
    '     Dim i, stepl, temp, pulse, number_steps_taken As Integer
    '     Dim ulstat As Integer
    '     ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)

    '     ulstat = cbAIn(BoardNum, 4, BIP10VOLTS, temp)

    '     For stepl = 1 To Math.Abs(number_of_steps_to_move)
    '         number_steps_taken = number_steps_taken + 1
    '         pulse = 1
    '         D_Out_Mask = &H80 * pulse + &H40 * direction
    '         ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
    '         For i = 1 To holdoff

    '         Next

    '         Hold_on_There(0.002)
    '         pulse = 0
    '         D_Out_Mask = &H80 * pulse + &H40 * direction
    '         ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)
    '         ulstat = cbToEngUnits(BoardNum, BIP10VOLTS, temp, limit_switch)
    '         If limit_switch < 0 Or number_steps_taken > 7000 Then
    '             Label13.Text = "out of range " & Str(limit_switch)
    '             geoffs_detector_position = geoffs_detector_position + (direction * number_steps_taken)
    '             Label13.Text = "pos. = " & Str(geoffs_detector_position)
    '             Exit Sub
    '         End If

    '     Next stepl
    '     geoffs_detector_position = geoffs_detector_position + (direction * number_steps_taken)
    '     Label13.Text = "pos. = " & Str(geoffs_detector_position)
    ' End Sub

    ' Private Sub Button15_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles movelin.Click
    '     mobile_detector(Int(mddir.Text)) ', Int(mdsteps.Text))
    ' End Sub


    ' Sub mobile_detector(ByVal move_cmd As Integer) ', ByVal holdoff As Integer)

    '     Dim ulstat As Integer
    '     move_cmd = Val(move_cmd)
    '     ulstat = cbDConfigPort(BoardNum, AUXPORT, DIGITALOUT)
    '     If move_cmd = 0 Then
    '         D_Out_Mask = (D_Out_Mask And 63) ' - &H40))
    '         '            D_Out_Mask = (D_Out_Mask And (&HFF - &H80))

    '     ElseIf move_cmd = 1 Then
    '         D_Out_Mask = (D_Out_Mask And 63)
    '         D_Out_Mask = (D_Out_Mask Or &H40)
    '     ElseIf move_cmd = 2 Then
    '         D_Out_Mask = (D_Out_Mask And 63)
    '         D_Out_Mask = (D_Out_Mask Or &H80)

    '     End If
    '     List1.Items.Add("mobile: " + Str(move_cmd) + " : " + Str(D_Out_Mask))

    '     ulstat = cbDOut(BoardNum, AUXPORT, D_Out_Mask)



    '     Hold_on_There(3)
    ' End Sub
    'Private Sub recalibrate_lambda_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles recalibrate_lambda.Click
    '    current_lambda = Val(LambdaText.Text)
    'My.Settings.lambda = current_lambda
    'My.Settings.Save()
    'End Sub

    
    Public Sub auto_gain_n(ByRef Current_Protocol, ByRef M_Number_Measuring_Lights, ByRef M_Number_Loops, ByRef M_Number_Pulses, _
                           ByRef Current_Intensity, ByRef m_Intensity, ByRef M_Take_Data, ByRef M_Primary_Gain, ByRef M_Reference_Gain, ByRef Measuring_Pulse_Duration, _
                           ByRef L_Measuring_Interval, ByRef Ref_Channel_Number, ByRef In_Channel, ByRef M_Measuring_Light, _
                           ByRef M_choke_actinic, ByVal Gain, ByVal points, _
                           ByRef data_time, ByRef Xe_Flash, ByRef q_switch, ByRef M_Far_Red, ByRef M_Blue_Actinic, ByRef Pre_Pulse, _
                           ByRef Pre_Pulse_Time, ByRef Pre_Delay, ByVal File_Name, ByRef exclude_from_gain)

        'Dim temp_current_protocol, temp_current_trace As Integer
        Dim primary_gain_index, use_primary_gain, i, ii, iii, M_Temp As Integer
        Dim zero_protocol As Integer = 0
        ' the integration times will increase as per the values in measuring_pulse_duration_arry
        ' Note: when used, the array is indexed to start with 1, so the first point is skipped

        Dim measuring_pulse_duration_array() As String = {"7u", "13u", "25u", "50u", "100u", "150u", "200u"} ' use only indexes 1-6
        Dim primary_gains() As Short = {0, 1, 2, 3} 'use index 1-4

        Dim measuring_pulse_duration_numerical_array(measuring_pulse_duration_array.Length() + 1) As Single
        Dim set_gain_for_this_measuring_light As Integer

        Gain(zero_protocol) = Gain(Current_Protocol) 'set the gain for the czero protocol to that in the current protocol

        Dim mi_index As Integer

        'go through each measuring light in the protocol being set.
        Dim best_times_array(M_Number_Measuring_Lights(Current_Protocol)) As String
        Dim best_gains_array(M_Number_Measuring_Lights(Current_Protocol)) As Integer
        Dim best_meas_gains_array(M_Number_Measuring_Lights(Current_Protocol)) As Integer
        Dim best_ref_gains_array(M_Number_Measuring_Lights(Current_Protocol)) As Integer

        For mi_index = 1 To M_Number_Measuring_Lights(Current_Protocol)
            Jabber_jabber.gain_data_list.Items.Add("------------------------------ Measuring light: " & Str(M_Measuring_Light(Current_Protocol, mi_index)) + " ------------------------------")
            If m_set_auto_gain(Current_Protocol, mi_index) > 0 Then 'only set the gains in the protocol if the defaul setting for m_set_auto_gain is 1

                set_gain_for_this_measuring_light = M_Measuring_Light(Current_Protocol, mi_index)

                For primary_gain_index = 0 To 3
                    use_primary_gain = primary_gains(primary_gain_index)
                    ' Dim number_primary_gains As Integer = detector_gain_array.Length() - 1
                    M_Number_Measuring_Lights(zero_protocol) = measuring_pulse_duration_array.Length() - 1 '* number_primary_gains
                    M_Number_Loops(zero_protocol) = 1
                    M_Number_Pulses(zero_protocol, 1) = 1 'fire each 'measuring light' only once.

                    'M_Intensity(zero_protocol, 0) = M_Intensity(Current_Protocol, 1)
                    m_Intensity(zero_protocol, 1) = Current_Intensity ' M_Intensity(Current_Protocol, 1)
                    M_Take_Data(zero_protocol, 1) = 1  ' want to take data, use default values

                    iii = 0 'usaed to calculate the index for the gains
                    For i = 1 To measuring_pulse_duration_array.Length() - 1
                        iii = iii + 1

                        M_Primary_Gain(zero_protocol, iii) = use_primary_gain
                        M_Reference_Gain(zero_protocol, iii) = use_primary_gain


                        parsenum(measuring_pulse_duration_array(i))
                        measuring_pulse_duration_numerical_array(i) = Time_Interval_in_Seconds
                        Measuring_Pulse_Duration(zero_protocol, iii) = measuring_pulse_duration_array(i)
                        'gain_data_list.Items.Add("mpd: " & Measuring_Pulse_Duration(zero_protocol, iii))

                        'Measuring_Pulse_Duration(zero_protocol, i) = measuring_pulse_duration_array(i)
                        L_Measuring_Interval(zero_protocol, iii, M_Number_Loops(zero_protocol)) = "1m"
                    Next

                    Current_Trace = 0
                    Ref_Channel_Number(zero_protocol) = Ref_Channel_Number(Current_Protocol)
                    In_Channel(zero_protocol) = In_Channel(Current_Protocol)

                    'the following sets in the zero protocol to the one being tested.
                    ' apparently there are 100 of these???

                    For M_Temp = 1 To M_Number_Measuring_Lights(zero_protocol)
                        M_Measuring_Light(zero_protocol, M_Temp) = set_gain_for_this_measuring_light
                    Next M_Temp


                    Call Multi_Trace(zero_protocol, M_Number_Loops, m_Intensity, M_Number_Pulses, M_Number_Measuring_Lights, _
                                     M_Measuring_Light, M_choke_actinic, L_Measuring_Interval, M_Primary_Gain, M_Reference_Gain, _
                                     Measuring_Pulse_Duration, Gain, In_Channel, _
                                     Ref_Channel_Number, points, _
                                     data_time, Xe_Flash, q_switch, M_Far_Red, M_Blue_Actinic, Pre_Pulse, Pre_Pulse_Time, Pre_Delay, _
                                     M_Take_Data)

                    Dim xxx As Integer = data_time.Length()

                    Auto_Gain_Calc(Current_Protocol, measuring_pulse_duration_numerical_array, points, File_Name(Current_Trace), _
                                   M_Number_Measuring_Lights(Current_Protocol), _
                                    In_Channel, Ref_Channel_Number, _
                                    measuring_pulse_duration_array.Length(), Measuring_Pulse_Duration, exclude_from_gain, _
                                    invert_raw, use_primary_gain)


                Next 'done cycling through the primary gains

                'Now we need find the gain settings that give the best (closest) to, but not over, the gain_slop setting
                Dim this_channel As Integer
                ' Importantly, we have to reconcile the two sets of gain settings because the same integration time is used for both channels.
                ' First, what is the highest integration time setting where there are valid (non-saturating) values for both channels
                ' for at least one primary gain
                ' For this_channel = 0 To 1 ' this_channel=0 means it is the samples channel data
                ' this_channel=0 means it is the samples channel data

                Dim highest_acceptable_integration_time_index = 0
                Dim highest_ok = 0

                'cycle through the primary gains
                ' this_channel=1 means it is the reference channel data
                'For primary_gain_index = 0 To 3
                '    any_not_ok = 0
                '    ' this_channel=1 means it is the reference channel data
                '    highest_acceptable_integration_time_index = M_Number_Measuring_Lights(zero_protocol) 'start with the maximum setting

                '    For M_Temp = 1 To M_Number_Measuring_Lights(zero_protocol) 'cycle through the integration times (which are stored as 
                '        Jabber_jabber.gain_data_list.Items.Add("primary_gain: " & Str(use_primary_gain) & " int: " & Measuring_Pulse_Duration(zero_protocol, M_Temp) & " Vs: " & Gain_Volts(0, use_primary_gain, M_Temp) & " Vr: " & Gain_Volts(1, use_primary_gain, M_Temp))
                '        If (Gain_Volts(0, use_primary_gain, M_Temp) < gain_slop) And (Gain_Volts(1, use_primary_gain, M_Temp) < gain_slop) Then
                '            highest_acceptable_integration_time_index = M_Temp 'if, at any integration times there are valid voltages, the set the flag to this index
                '        End If
                '    Next
                'Next


                ' this_channel=1 means it is the reference channel data
                highest_acceptable_integration_time_index = M_Number_Measuring_Lights(zero_protocol) 'start with the maximum setting
                Dim chan0_best_primary_gain As Integer = 0
                Dim chan1_best_primary_gain As Integer = 0

                Dim chan0_ok, chan1_ok As Integer

                For M_Temp = 1 To M_Number_Measuring_Lights(zero_protocol)  'cycle through the integration times starting from the lowest
                    highest_ok = -1 'start at negative 1 
                    chan0_ok = -1
                    chan1_ok = -1

                    For primary_gain_index = 0 To 3 ' go through all the primary gains from low to high 
                        Jabber_jabber.gain_data_list.Items.Add("primary_gain: " & Str(primary_gain_index) & " int: " & Measuring_Pulse_Duration(zero_protocol, M_Temp) & " Vs: " & Gain_Volts(0, primary_gain_index, M_Temp) & " Vr: " & Gain_Volts(1, primary_gain_index, M_Temp))
                        If (Gain_Volts(0, primary_gain_index, M_Temp) < gain_slop) Then ' if the value is less than the max, save the gain index in chan0_ok
                            chan0_ok = primary_gain_index
                        End If

                        If (Gain_Volts(1, primary_gain_index, M_Temp) < gain_slop) Then ' if the value is less than the max, save the gain index in chan1_ok
                            chan1_ok = primary_gain_index
                            'highest_ok = primary_gain_index 'if something works, 
                            'highest_acceptable_integration_time_index = M_Temp 'if, at any integration times there are valid voltages, the set the flag to this index
                        End If
                    Next

                    If (chan0_ok > -1 And chan1_ok > -1) Then 'if there are good values for each channel then save the primary gains for each channel in chan0_best and chan1_best
                        chan0_best_primary_gain = chan0_ok
                        chan1_best_primary_gain = chan1_ok
                        highest_acceptable_integration_time_index = M_Temp 'if, at any integration times there are valid voltages, the set the flag to this index
                    End If


                Next



                Jabber_jabber.gain_data_list.Items.Add("highest OK int time" & Str(this_channel) & " is " & Measuring_Pulse_Duration(zero_protocol, highest_acceptable_integration_time_index))
                ' now, for each channel, find the highest acceptable primary gain 

                Measuring_Pulse_Duration(Current_Protocol, mi_index) = Measuring_Pulse_Duration(zero_protocol, highest_acceptable_integration_time_index)

                Jabber_jabber.gain_data_list.Items.Add("highest OK sample primary gain is :" & Str(chan0_best_primary_gain))
                Jabber_jabber.gain_data_list.Items.Add("highest OK reference primary gain is :" & Str(chan1_best_primary_gain))

                M_Primary_Gain(Current_Protocol, mi_index) = chan0_best_primary_gain
                M_Reference_Gain(Current_Protocol, mi_index) = chan1_best_primary_gain


                '    If this_channel = 0 Then
                '        M_Primary_Gain(Current_Protocol, mi_index) = highest_acceptable_primary_gain(this_channel)

                '    Else

                '    End If


                '    Dim highest_acceptable_primary_gain() = {0, 0}
                '    'Dim ttest As String
                '    Dim this_channel_number As Short
                '    Dim this_channel_name As String
                '    For this_channel = 0 To 1
                '        If this_channel = 0 Then
                '            this_channel_name = "sample"
                '            this_channel_number = In_Channel(Current_Protocol)
                '        Else
                '            this_channel_name = "reference"
                '            this_channel_number = Ref_Channel_Number(Current_Protocol)
                '        End If
                '        For primary_gain_index = 0 To 3 'cycle through the primary gains
                '            use_primary_gain = primary_gains(primary_gain_index)
                '            'ttest = "LED# :" + Str(mi_index) + " " + this_channel_name + ": " + Str(this_channel_number) + ":" + Str(primary_gain_index)
                '            'Jabber_jabber.gain_data_list.Items.Add(ttest)

                '            'For use_primary_gain = 0 To ok_primary_gains.Length() 'cycle through the primary gains
                '            If (Gain_Volts(this_channel, use_primary_gain, highest_acceptable_integration_time_index) < gain_slop) Then
                '                highest_acceptable_primary_gain(this_channel) = use_primary_gain
                '            End If
                '        Next
                '        Jabber_jabber.gain_data_list.Items.Add("high OK primary gain channel " & Str(this_channel) & " is " & Str(highest_acceptable_primary_gain(this_channel)))
                '        If this_channel = 0 Then
                '            M_Primary_Gain(Current_Protocol, mi_index) = highest_acceptable_primary_gain(this_channel)

                '        Else
                '            M_Reference_Gain(Current_Protocol, mi_index) = highest_acceptable_primary_gain(this_channel)

                '        End If
                '    Next



            End If

        Next

    End Sub

    Private Sub zg5_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles zg5.Load

    End Sub

    Private Sub zg4_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles zg4.Load

    End Sub

    Private Sub MonochromatorToolStripMenuItem_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MonochromatorToolStripMenuItem.Click
        Form3.Show()
    End Sub

    Private Sub Button1_Click_3(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button1.Click
        Form3.moveTo.Text = TextBox4.Text
        Label14.Text = "moving"

        Form3.moveIt()

        While Form3.moving = True
            System.Windows.Forms.Application.DoEvents()
        End While

        Label14.Text = "done"

    End Sub
End Class