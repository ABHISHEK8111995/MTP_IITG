function varargout = non_linear_vib(varargin)
% NON_LINEAR_VIB MATLAB code for non_linear_vib.fig
%      NON_LINEAR_VIB, by itself, creates a new NON_LINEAR_VIB or raises the existing
%      singleton*.
%
%      H = NON_LINEAR_VIB returns the handle to a new NON_LINEAR_VIB or the handle to
%      the existing singleton*.
%
%      NON_LINEAR_VIB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NON_LINEAR_VIB.M with the given input arguments.
%
%      NON_LINEAR_VIB('Property','Value',...) creates a new NON_LINEAR_VIB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before non_linear_vib_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to non_linear_vib_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help non_linear_vib

% Last Modified by GUIDE v2.5 08-Aug-2023 01:02:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @non_linear_vib_OpeningFcn, ...
                   'gui_OutputFcn',  @non_linear_vib_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before non_linear_vib is made visible.
function non_linear_vib_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to non_linear_vib (see VARARGIN)

% Choose default command line output for non_linear_vib
% handles.output = hObject;
handles.output = hObject;
Unt=get(0,'Units');
p=get(0,'Screensize');
set(hObject,'Units',Unt);
figsize=get(hObject,'Position');
wf=figsize(3);
hf=figsize(4);ws=p(3);
hs=p(4);
pfig=[(ws-wf)/2 (hs-hf)/2 wf hf];
set(hObject,'position',pfig);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes non_linear_vib wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = non_linear_vib_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global den
den=str2double(get(handles.edit1,'string'));




% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global d
d=str2double(get(handles.edit2,'string'));



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

global l
l=str2double(get(handles.edit3,'string'));


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

global g
g=str2double(get(handles.edit4,'string'));


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

global M
M=str2double(get(handles.edit5,'string'));


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

global R
R=str2double(get(handles.edit6,'string'));


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

global alpha
alpha=str2double(get(handles.edit7,'string'));


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

global f
f=str2double(get(handles.edit8,'string'));



% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

global zeta
zeta=str2double(get(handles.edit9,'string'));



% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

global x y
x=(get(handles.radiobutton1,'value'));
y=set(handles.radiobutton2,'value',0);
y=0;

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

global y x
y=(get(handles.radiobutton2,'value'));
x=set(handles.radiobutton1,'value',0);
x=0;


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

global amin
amin=str2double(get(handles.edit13,'string'));


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

global amax
amax=str2double(get(handles.edit14,'string'));



% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x y amin amax den d l i_len alpha f zeta g M R f_len freq_diff1

format long g

vol=(pi/4)*d^2*l;
m=vol*den;
i=M*R^2/2;
D0=[];   %this is to make it empty if there is no optimise dia so that errordlg works
if x==1
    
%     if (~isfinite(freq_diff)| isempty(freq_diff))
%         errordlg('PLEASE PROVIDE DIFFERENCE IN FREQUENCY IN CONTROL SETTING','ERROR','modal')
%         return
%     else 
%     end
      
    k=0.99;  %initial guess
    di=inf;  
    j=0;
    r=0;
    j_s=pi*d^4/32;
    i_s=(i+m*d^2/24);
    f_s=sqrt(g*j_s/(l*i_s));
    
   while(k>0 & di>0.1*d)   %di is decided by manufacturing possiblity of making hole in shaft
    q=sqrt(1+1/k^2);  %q=d0/di
    di=k*d;
    d0=q*di;
    j_h=pi*(d0^4-di^4)/32;
    i_h=(i+m*(d0^2+di^2)/24);
    j=j+1;  %for no.of iterations
    f_h=sqrt(g*j_h/(l*i_h));
    F=abs(sqrt(j_s/i_s)-sqrt(j_h/i_h))/sqrt(j_s/i_s);

    if (F*100<=freq_diff1)     
        k=k-0.005;
        r=r+1;
        D0(r)=d0;
        Di(r)=di;
        freq1(r)=f_h;            
    else
        k=k-0.005;
        
    end
    
end 
    
    if isempty(D0)
        warndlg('THERE IS NO OPTIMISE DIAMETER FOR THE GIVEN INPUT DATA','WARNING','modal') 
        return
    else
    end
    
         D=[D0;Di;freq1]';
         p=size(D,1);
            for i1=1:p
                q(i1)=i1;
            end 
            
        a=[q',D];          
        [s,v] = listdlg('PromptString','Select the figures','Name','Figures Options',...
                'SelectionMode','multiple','ListSize',[250 300],...
                'ListString',num2str(a));
            
            if v==1
                  j_h1=pi*(D0(s)^4-Di(s)^4)/32;  
                  i_h1=(i+m*(D0(s)^2+Di(s)^2)/24);     
            else
                return
            end
     f_h=sqrt(g*j_h1/(l*i_h1));
    
elseif y==1
    
    i_len1=i_len;
    k=0.99;  %initial guess
    r1=0;               
    F1=modal_freq_solid(d,l,g,den,i);
    f_s=F1;
    while (i_len1<f_len)   %variable    %di is decided by manufacturing possiblity of making hole in shaft
                                       %now it is 10 % of initial dia and initial length 
        while(k>0)                                      
            di=k*d;
            D=sqrt(d^2+(k^2*d^2*l)/i_len1);      %d0=d
             F2=modal_freq_step(d,di,D,l,i_len1,g,den,i);
             F=abs(F1-F2)/F1;
            if (F*100<=freq_diff1)
                k=k-0.005;            %variable
                r1=r1+1;
                D0(r1)=D;
                Di(r1)=di;
                i_length(r1)=i_len1;
                freq1(r1)=F2;       %angular frequency
            else
                k=k-0.005;
            end
        end

        i_len1=i_len1+0.05*l;
        k=0.99;
    % dia_original=sqrt(((d^2-di^2)*(l-i_len)+(D^2-di^2)*i_len)/l); %this is mass conservation  
    end 
    
     if isempty(D0)
        warndlg('THERE IS NO OPTIMISE DIAMETER FOR THE GIVEN INPUT DATA','WARNING','modal') 
        return
    else
    end

    D=[D0;Di;i_length;freq1]';     % optimise_diameter
    
    p=size(D,1);
         q=[];
            for i1=1:p
                q(i1)=i1;
            end 
            
        a=[q',D];          
        [s,v] = listdlg('PromptString','Select the figures','Name','Figures Options',...
                'SelectionMode','multiple','ListSize',[250 300],...
                'ListString',num2str(a));
            
            if v==1
                  f_h=freq1(s);    
            else
                return
            end
    
else
   errordlg('FIRST SELECT THE PRISMATIC/NON PRISMATIC BUTTON','ERROR','modal');
   return
end

    a=amin:0.01:amax;
    sigma_s1=(3/8)*(alpha*a.^2/f_s) + sqrt((f^2./(4*f_s^2*a.^2))-zeta^2);
    sigma_s2=(3/8)*(alpha*a.^2/f_s) - sqrt((f^2./(4*f_s^2*a.^2))-zeta^2);
    sigma_h1=(3/8)*(alpha*a.^2/f_h) + sqrt((f^2./(4*f_h^2*a.^2))-zeta^2);
    sigma_h2=(3/8)*(alpha*a.^2/f_h) - sqrt((f^2./(4*f_h^2*a.^2))-zeta^2);
    
    figure()
    plot(sigma_s1,a,'--g',sigma_s2,a,'--r','Linewidth',2)
    hold on
    plot(sigma_h1,a,'b',sigma_h2,a,'k','Linewidth',1)
    grid minor 
    xlabel('SIGMA')
    ylabel('AMPLITUDE')
    title('FREQUENCY RESPONSE OF NON LINEAR ROTOR SHAFT SYSTEM')
    legend('Original shaft','Original shaft','Stepped Shaft','Stepped Shaft')

 
    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
clear all


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

global i_len
i_len=str2double(get(handles.edit11,'string'));


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

global f_len
f_len=str2double(get(handles.edit12,'string'));



% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over figure background.




function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double

global freq_diff1
freq_diff1=str2double(get(handles.edit15,'string'));




% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
