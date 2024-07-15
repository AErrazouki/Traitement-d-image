function varargout = Tp(varargin)
% TP MATLAB code for Tp.fig
%      TP, by itself, creates a new TP or raises the existing
%      singleton*.
%
%      H = TP returns the handle to a new TP or the handle to
%      the existing singleton*.
%
%      TP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TP.M with the given input arguments.
%
%      TP('Property','Value',...) creates a new TP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tp

% Last Modified by GUIDE v2.5 24-Mar-2024 20:18:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tp_OpeningFcn, ...
                   'gui_OutputFcn',  @Tp_OutputFcn, ...
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


% --- Executes just before Tp is made visible.
function Tp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tp (see VARARGIN)

% Choose default command line output for Tp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Fichier_Callback(hObject, eventdata, handles)
% hObject    handle to Fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Bruit_Callback(hObject, eventdata, handles)
% hObject    handle to Bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Transformation_Callback(hObject, eventdata, handles)
% hObject    handle to Transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Filtre_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Filtre_Freq_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre_Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to Ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.ImgO)
handles.courant_data = handles.ima;
subimage(handles.courant_data);
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Quiter_Callback(hObject, eventdata, handles)
% hObject    handle to Quiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --------------------------------------------------------------------
function Enregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to Enregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');


% --------------------------------------------------------------------
function Gauss_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'gaussian');
axes(handles.ImgT)
imshow(Imgb);
handles.courant_data=Imgb;
handles.output=hObject;

% --------------------------------------------------------------------
function Poivre_Sel_Callback(hObject, eventdata, handles)
% hObject    handle to Poivre_Sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'Salt & Pepper', 0.01);
axes(handles.ImgT)
imshow(Imgb);
handles.courant_data=Imgb;
handles.output=hObject;


% --------------------------------------------------------------------
function Inv_coul_Callback(hObject, eventdata, handles)
% hObject    handle to Inv_coul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image = double(image);
[l c]=size(image);
v=image;
for i=1:l
   for j=1:c
     v(i,j)=255-double(image(i,j));
    end
 end 
v=uint8(v); 
axes(handles.ImgT);
subimage(v);

% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to Histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
d = length(size(img));
if d==3
    I = rgb2gray(img);
elseif d==2
     I = img
end
axes(handles.ImgO);
subimage(I);
%[H n]=imhist(I)
[nl nc]=size(I);
v=double(I);
vec=[1:256];
l=0;
for k=0:255 
   for i=1:nl
        for j=1:nc
            if v(i,j)==k 
              l=l+1;
           end
        end
    end
   vec(k+1)=l;
   l=0;
end
axes(handles.ImgT);
plot(vec)
%plot(H);


% --------------------------------------------------------------------
function Contraste_Callback(hObject, eventdata, handles)
% hObject    handle to Contraste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Lecture de l'image à partir de handles.courant_data
image = handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
  
image = double(image);
A = double(min(min(image)));
B = double(max(max(image)));
P = 255 / (B - A);
L = -P * A;
[l, c] = size(image);
v = image;
for i = 1:l
    for j = 1:c
        fpixel = image(i,j) * P + L;
         if(fpixel > 255)
             fpixel = 255;
         elseif(fpixel < 0)
             fpixel = 0;
         end
          v(i,j) = fpixel;
     end
end
v = uint8(v);
axes(handles.ImgT);
imshow(v);
figure, imhist(v,256);


% --------------------------------------------------------------------
function Luminosite_Callback(hObject, eventdata, handles)
% hObject    handle to Luminosite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
    for j=1:c
        pix=image(i,j)+50;
        if(pix>255)
           pix=255;
        elseif (pix<0)
            pix=0;
        end
        v(i,j)=pix;    
    end
end  
v=uint8(v); 
axes(handles.ImgT);
subimage(v);

% --------------------------------------------------------------------
function Binarisation_Callback(hObject, eventdata, handles)
% hObject    handle to Binarisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I4 = handles.courant_data;
if ndims(I4) == 3
    I4 = rgb2gray(I4);
end
image = double(I4);
[l c]=size(image);
seuil = 128;
for i=1:l
    for j=1:c
            if image(i, j) > seuil
                image(i, j) = 255;
            else
               image(i, j) = 0;
            end
    end
end
axes(handles.ImgT);
imshow(image);


% --------------------------------------------------------------------
function Niv_gris_Callback(hObject, eventdata, handles)
% hObject    handle to Niv_gris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ima=handles.courant_data;
d = ndims(ima);
if d==3
    imagray=rgb2gray(ima);
elseif d==2
   imagray=ima;
end
axes(handles.ImgT);
subimage(imagray);


% --------------------------------------------------------------------
function Fpb_ideal_Callback(hObject, eventdata, handles)
% hObject    handle to Fpb_ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data;
F=fftshift(fft2(I));  
[m, n] = size(F);
D0 = 10;
H = zeros(m, n);
for u = 1:m
    for v = 1:n
        D = sqrt((u - m/2)^2 + (v - n/2)^2);
        if D <= D0
            H(u, v) = 1;
        end
    end
end
for i=1:m 
   for j=1:n 
     G(i,j)=F(i,j)*H(i,j); 
   end 
end 

g=ifft2(G);
axes(handles.ImgT);
imshow(abs(g),[0,255]);

% --------------------------------------------------------------------
function Lineaire_Callback(hObject, eventdata, handles)
% hObject    handle to Lineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Non_Lineaire_Callback(hObject, eventdata, handles)
% hObject    handle to Non_Lineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Moyenneur_Callback(hObject, eventdata, handles)
% hObject    handle to Moyenneur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gaussien_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gauss_3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss_3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/16)*[1 2 1 ;2 4 2 ; 1 2 1];
for x = 2 : n-1
    for y = 2 : m-1
        f=image(x-1:x+1,y-1:y+1);
        v=f.*H;
        b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Gauss_5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Gauss_5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/256)*[1 4 6 4 1 ; 4 16 24 16 4 ; 6 24 36 24 6 ; 4 16 24 16 4 ; 1 4 6 4 1];
for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
       v=f.*H;
       b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Moy_3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Moy_3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/9)*[1 1 1 ; 1 1 1 ; 1 1 1 ];
for x = 2 : n-1
    for y = 2 : m-1
      f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);


% --------------------------------------------------------------------
function Moy_5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Moy_5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1];
for x = 3 : n-2
    for y = 3 : m-2
      f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);


% --------------------------------------------------------------------
function Mediane_Callback(hObject, eventdata, handles)
% hObject    handle to Mediane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image=double(image);
[n,m]=size(image);
img=image;
for i=2:n-1
    for j=2:m-1
       fenetre=image(i-1:i+1,j-1:j+1);
       v=[fenetre(1,:) fenetre(2,:) fenetre(3,:)];
       v=sort(v);
       a=median(v);
       img(i,j)=a;
    end
end
b=uint8(img);
handles.ima_traite = b;
axes(handles.ImgT);
subimage(b);
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Conique_Callback(hObject, ~, handles)
% hObject    handle to Conique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[0 0 1 0 0 ; 0 2 2 2 0 ; 1 2 5 2 1 ; 0 2 2 2 0 ; 0 0 1 0 0];
for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
       v=f.*H;
       b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
%guidata(hObject, handles);

% --------------------------------------------------------------------
function Pyramidal_Callback(hObject, eventdata, handles)
% hObject    handle to Pyramidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/81)*[1 2 3 2 1 ; 2 4 6 4 2 ; 3 6 9 6 3 ; 2 4 6 4 2 ; 1 2 3 2 1];
for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
       v=f.*H;
       b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Fph_ideal_Callback(hObject, eventdata, handles)
% hObject    handle to Fph_ideal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data; 
F=fftshift(fft2(I)); 
[m, n] = size(F);
D0 = 1;
H = zeros(m, n);
for u = 1:m
    for v = 1:n
        D = sqrt((u - m/2)^2 + (v - n/2)^2);
        if D >= D0
            H(u, v) = 1;
        end
    end
end
G = F .* H;
g=ifft2(G);
axes(handles.ImgT);
imshow(255-abs(g),[0,255]);


% --------------------------------------------------------------------
function Laplcien_Callback(hObject, eventdata, handles)
% hObject    handle to Laplcien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
[n,m]=size(image);
image = double(image);
b=zeros(n,m);
M1=[-1 -1 -1;-1 8 -1;-1 -1 -1];
for i=2:n-1
    for j=2:m-1
        V=image((i-1:i+1),(j-1:j+1));
        S=V.*M1;
        b(i,j)=sum(S(:));
    end
end
b=uint8(b);
axes(handles.ImgT);
subimage(b);


% --------------------------------------------------------------------
function Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
[n,m]=size(image);
image = double(image);
output=image;
[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [0,0,0;-1,0,1;0,0,0]; 
maskver = [0,-1,0;0,0,0;0,1,0];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end
for i=4:(m-3)
    for j=4:(n-3)       
       output(i,j)=sqrt(outputhor(i,j).^2 + outputver(i,j)^2);
    end 
end 
output=uint8(output); 
axes(handles.ImgT);
subimage(output);


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
[n,m]=size(image);
image = double(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [-1,0,1;-2,0,2;-1,0,1]; 
maskver = [-1,-2,-1;0,0,0;1,2,1];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
for i=4:(m-3)
   for j=4:(n-3) 
        output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
 output=uint8(output); 
axes(handles.ImgT);
subimage(output);



% --------------------------------------------------------------------
function Prewit_Callback(hObject, eventdata, handles)
% hObject    handle to Prewit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
[n,m]=size(image);
image = double(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [-1,0,1;-1,0,1;-1,0,1]; 
maskver = [-1,-1,-1;0,0,0;1,1,1];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
end
for i=4:(m-3)
   for j=4:(n-3) 
       output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
for i=1:n-1
    for j=1:m-1
        if output(i,j) < 30
            output(i,j) =0;
        end
     end
end
output=uint8(output); 
axes(handles.ImgT);
subimage(output);


% --------------------------------------------------------------------
function Robert_Callback(hObject, eventdata, handles)
% hObject    handle to Robert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
[n,m]=size(image);
image = double(image);
for x=1:n-1
  for y=1:m-1
    b(x,y)= abs(uint8( double(image(x,y))-double(image(x+1,y+1))))+ abs(uint8( double(image(x,y+1)) - double(image(x+1,y))));
  end
end
%Seuillage
for i=1:n-1
    for j=1:m-1
        if b(i,j) < 30
            b(i,j)=0;
        end
     end
end
handles.ima_traite = b;
axes(handles.ImgT);
subimage(b);


% --------------------------------------------------------------------
function Fpb_butter_Callback(hObject, eventdata, handles)
% hObject    handle to Fpb_butter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data;
F=fftshift(fft2(I));    
[M,N]=size(F);
H0=zeros(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H0(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
n=3; 
for i=1:M 
    for j=1:N  
        H(i,j)=1/(1+(H0(i,j)/D0)^(2*n)); 
        G(i,j)=F(i,j)*H(i,j); 
    end 
end 
g=ifft2(G); 
axes(handles.ImgT);
imshow(abs(g),[0,255]);

% --------------------------------------------------------------------
function Fph__butter_Callback(hObject, eventdata, handles)
% hObject    handle to Fph__butter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
F=fftshift(fft2(I));  
[M,N]=size(F); 
H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
n=3; 
for i=1:M 
    for j=1:N 
        H(i,j)=1/(1+(D0/H1(i,j))^(2*n)); 
        G(i,j)=F(i,j)*H(i,j);  
    end 
end 
g=ifft2(G); 
axes(handles.ImgT);
imshow(255-abs(g),[0,255]);


% --------------------------------------------------------------------
function Coin_Callback(hObject, eventdata, handles)
% hObject    handle to Coin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Susan_Callback(hObject, eventdata, handles)
% hObject    handle to Susan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
d = length(size(image));
if d == 3
    image = double(rgb2gray(image));
elseif d == 2
    image = double(image);
end
[n, m] = size(image);
rayon = 2;
alpha = 50;
r = 2;
alpha = alpha / 100;
% Générateur de masque
mask = zeros(2 * rayon + 1);
b = ones(rayon + 1);
for i = 1:rayon + 1
   for j = 1:rayon + 1
        if (rayon == 1)
             if (j > i)
                 b(i, j) = 0;
             end
         else
             if (j > i + 1)
                 b(i, j) = 0;
             end
         end
    end
end
mask(1:rayon + 1, rayon + 1:2 * rayon + 1) = b;
mask(1:rayon + 1, 1:rayon + 1) = rot90(b);
mask0 = mask;
mask0 = flipdim(mask0, 1);
mask = mask0 + mask;
mask(rayon + 1, :) = mask(rayon + 1, :) - 1;
% Réponse maximale
max_response = sum(mask(:));
% Balayage de toute l'image
 f = zeros(n, m);
    
    for i = (rayon + 1):n - rayon
        for j = (rayon + 1):m - rayon
            image_courant = image(i - rayon:i + rayon, j - rayon:j + rayon);
            image_courant_mask = image_courant .* mask;
            intensity_central = image_courant_mask(rayon + 1, rayon + 1);
            s = exp(-1 * (((image_courant_mask - intensity_central) / max_response).^6));
            somme = sum(s(:));
            
            % Soustraire les zéros des filtres si le centre du masque est un 0
            if (intensity_central == 0)
                somme = somme - length(find(mask == 0));
            end
            
            f(i, j) = somme;
        end
    end
    
    % Sélection et seuillage des points d'intérêt
    ff = f(rayon + 1:n - (rayon + 1), rayon + 1:m - (rayon + 1));
    minf = min(min(ff));
    maxf = max(max(f));
    fff = f;
    d = 2 * r + 1;
    temp1 = round(n / d);
    
    if (temp1 - n / d) < 0.5 && (temp1 - n / d) > 0
        temp1 = temp1 - 1;
    end
    
    temp2 = round(m / d);
    
    if (temp2 - m / d) < 0.5 && (temp2 - m / d) > 0
        temp2 = temp2 - 1;
    end
    
    fff(n:temp1 * d + d, m:temp2 * d + d) = 0;
    seuil = minf + alpha * (maxf - minf);
    [u, v] = find(minf <= fff & fff <= seuil);
    
    % Affichage des résultats
    axes(handles.ImgT);
    imshow(image);
    hold on;
    plot(v, u, '.r', 'MarkerSize', 10);
    hold off;

% --------------------------------------------------------------------
function Harris_Callback(hObject, eventdata, handles)
% hObject    handle to Harris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
%==========================================================================
if(size(img,3)==3)
    img=rgb2gray(img);
end
%==========================================================================
k=0.04;
sigma=1; seuil=200; r=6; w=5*sigma;
[m,n]=size(img)
imd=double(img);
dx=[-1 0 1 , -2 0 2 ,-1 0 1]; % derivée horizontale : filtre de Sobel
dy=[-1 -2 -1 , 0 0 0 ,1 2 1]; % derivée verticale : filtre de Sobel

g = fspecial('gaussian',max(1,fix(w)), sigma);

Ix=conv2(imd,dx,'same'); % Dérivée horizontale
Iy=conv2(imd,dy,'same');

Ix2=conv2(Ix.^2, g, 'same');
Iy2=conv2(Iy.^2, g, 'same');
Ixy=conv2(Ix.*Iy, g,'same');

detM=Ix2.*Iy2-Ixy.^2;
trM=Ix2+Iy2;
R=detM-k*trM.^2;
%==========================================================================
R1=(1000/(1+max(max(R))))*R;
%==========================================================================          
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
    R1(u(k),v(k))=0;
end

R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;

% Non-maximal suppression pour ne garder que les coins les plus saillants
[m1,n1]=size(R11);
for i=r+1:m1-r
    for j=r+1:n1-r
        fenetre=R11(i-r:i+r,j-r:j+r);
        ma=max(max(fenetre));
        if fenetre(r+1,r+1)<ma
            R11(i,j)=0;
        end
    end
end

R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);

nb = length(x);
axes(handles.ImgT);
imshow(img);
hold on;
plot(y, x, 'r.');
%title(['Nombre de coins détectés : ', num2str(nb)]);


% --------------------------------------------------------------------
function Morphologie_Callback(hObject, eventdata, handles)
% hObject    handle to Morphologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Dilatation_Callback(hObject, eventdata, handles)
% hObject    handle to Dilatation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
b=[1 1 1;1 1 1;1 1 1];
if size(image, 3) == 3
    image = rgb2gray(image);
end
a1=imdilate(image,b);
axes(handles.ImgT);
subimage(a1);


% --------------------------------------------------------------------
function Erosion_Callback(hObject, eventdata, handles)
% hObject    handle to Erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Convertir l'image en niveaux de gris si elle est en couleur
image = handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
b=[1 1 1;1 1 1;1 1 1];
a2=imerode(image,b);
axes(handles.ImgT);
subimage(a2);


% --------------------------------------------------------------------
function Ouverture_Callback(hObject, eventdata, handles)
% hObject    handle to Ouverture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
b=[1 1 1;1 1 1;1 1 1];
%a2=imerode(image,b);
%a1=imdilate(a2,b);
a1=imopen(image,b);
axes(handles.ImgT);
subimage(a1);


% --------------------------------------------------------------------
function Fermeture_Callback(hObject, eventdata, handles)
% hObject    handle to Fermeture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
    image = rgb2gray(image);
end
b=[1 1 1;1 1 1;1 1 1];
%a1=imdilate(image,b);
%a2=imerode(a1,b);
a2=imclose(image,b);
axes(handles.ImgT);
subimage(a2);


% --------------------------------------------------------------------
function Hough_Callback(hObject, eventdata, handles)
% hObject    handle to Hough (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Egalisation_Callback(hObject, eventdata, handles)
% hObject    handle to Egalisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
if ndims(img) == 3
    img = rgb2gray(img);
end

h = imhist(img);
n = h / sum(h);
cdf = cumsum(n);
img_eq = uint8(255 * cdf(img + 1));
axes(handles.ImgT);
imshow(img_eq)


% --------------------------------------------------------------------
function FPBande_Callback(hObject, eventdata, handles)
% hObject    handle to FPBande (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
     image = rgb2gray(image);
end

F = fftshift(fft2(image));
[M, N] = size(F);
% Définir les paramètres du filtre passe-bande
D0 = 50; % Fréquence de coupure basse
D1 = 150; % Fréquence de coupure haute
% Créer un masque pour le filtre passe-bande
H = zeros(M, N);
for u = 1:M
     for v = 1:N
        D = sqrt((u - M/2)^2 + (v - N/2)^2);
        if D >= D0 && D <= D1
             H(u, v) = 1;
        end
     end
end
G = F .* H;
g = ifft2(ifftshift(G));
axes(handles.ImgT);
imshow(abs(g), []);


% --------------------------------------------------------------------
function Gradient_M_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gradient_interne_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_interne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
if size(image, 3) == 3
     image = rgb2gray(image);
end
se = strel('disk', 3);
fond = imerode(image, se);
gradient_interne = image - fond;
axes(handles.ImgT);
imshow(gradient_interne);


% --------------------------------------------------------------------
function Gradient_externe_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_externe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
if size(image, 3) == 3
     image = rgb2gray(image);
end
se = strel('disk', 3);
dilated_image = imdilate(image, se);
gradient_externe = dilated_image - image;
axes(handles.ImgT);
imshow(gradient_externe);


% --------------------------------------------------------------------
function Gradient_morphologique_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_morphologique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
if size(image, 3) == 3
     image = rgb2gray(image);
end
se = strel('disk', 3);
dilated_image = imdilate(image, se);
eroded_image = imerode(image, se);
gradient_morphologique = dilated_image - eroded_image;
axes(handles.ImgT);
imshow(gradient_morphologique);


% --------------------------------------------------------------------
function Kirch_Callback(hObject, eventdata, handles)
% hObject    handle to Kirch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Récupérer les données de l'image à partir des handles
image = handles.courant_data;
x = double(image);
    
% Définir les masques de Kirsch pour les 8 directions (45°)
kirsch_masks = {
     [5,5,5; -3,0,-3; -3,-3,-3],
     [5,5,-3; 5,0,-3; -3,-3,-3],
     [5,-3,-3; 5,0,-3; 5,-3,-3],
     [-3,-3,-3; 5,0,-3; 5,5,-3],
     [-3,-3,-3; -3,0,-3; 5,5,5],
     [-3,-3,-3; -3,0,5; -3,5,5],
     [-3,-3,5; -3,0,5; -3,-3,5],
     [-3,5,5; -3,0,5; -3,-3,-3]
};
    
outputImage = zeros(size(image));

for i = 1:8
    conv_result = imfilter(x, kirsch_masks{i}, 'replicate');
    outputImage = max(outputImage, conv_result);
end
axes(handles.ImgT);
subimage(uint8(outputImage));


% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
if size(img, 3) == 3
    img = rgb2gray(img);
end

T_Low = 0.075;
T_High = 0.175;

B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;

A=conv2(img, B, 'same');

KGx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];

Filtered_X = conv2(A, KGx, 'same');
Filtered_Y = conv2(A, KGy, 'same');
%calcule des angles pour chaque pixel dans la direction dans laquelle la
%luminosite de l'image change plus rapidement
arah = atan2 (Filtered_Y, Filtered_X);
arah = arah*180/pi;

pan=size(A,1);
leb=size(A,2);

%Adjustment for negative directions, making all directions positive
for i=1:pan
    for j=1:leb
        if (arah(i,j)<0) 
            arah(i,j)=360+arah(i,j);
        end;
    end;
end;
arah2=zeros(pan, leb);

for i = 1  : pan
    for j = 1 : leb
        if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) || (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) || (arah(i, j) >= 337.5) && (arah(i, j) <= 360))
            arah2(i, j) = 0;
        elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) || (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
            arah2(i, j) = 45;
        elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) || (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
            arah2(i, j) = 90;
        elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) || (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
            arah2(i, j) = 135;
        end;
    end;
end;
%Calculate de changement d'intensite 
magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
magnitude2 = sqrt(magnitude);

BW = zeros (pan, leb);

%Non-Maximum Supression
for i=2:pan-1
    for j=2:leb-1
        %en fct de la dicrection de changement 
        if (arah2(i,j)==0)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
        elseif (arah2(i,j)==45)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
        elseif (arah2(i,j)==90)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
        elseif (arah2(i,j)==135)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
        end;
    end;
end;

BW = BW.*magnitude2;
%Hysteresis Thresholding
T_Low = T_Low * max(max(BW));
T_High = T_High * max(max(BW));

T_res = zeros (pan, leb);

for i = 1  : pan
    for j = 1 : leb
        if (BW(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (BW(i, j) > T_High)
            T_res(i, j) = 1;
        %Using 8-connected components
        elseif ( BW(i+1,j)>T_High || BW(i-1,j)>T_High || BW(i,j+1)>T_High || BW(i,j-1)>T_High || BW(i-1, j-1)>T_High || BW(i-1, j+1)>T_High || BW(i+1, j+1)>T_High || BW(i+1, j-1)>T_High)
            T_res(i,j) = 1;
        end;
    end;
end;

edge_final = uint8(T_res.*255);
axes(handles.ImgT);
subimage(edge_final);


% --------------------------------------------------------------------
function MarrHildreth_Callback(hObject, eventdata, handles)
% hObject    handle to MarrHildreth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
im=double(img);
gfilter= [0 0 1 0 0;
       0 1 2 1 0;
       1 2 -16 2 1;
       0 1 2 1 0;
       0 0 1 0 0];
   
smim=conv2(im,gfilter);
[rr,cc]=size(smim);
zc=zeros([rr,cc]);
for i = 2:rr-1
    for j = 2:cc-1
        if smim(i, j) > 0
            if (smim(i, j+1) >= 0 && smim(i, j-1) < 0) || (smim(i, j+1) < 0 && smim(i, j-1) >= 0)
                zc(i, j) = smim(i, j);
            elseif (smim(i+1, j) >= 0 && smim(i-1, j) < 0) || (smim(i+1, j) < 0 && smim(i-1, j) >= 0)
                zc(i, j) = smim(i, j);
            elseif (smim(i+1, j+1) >= 0 && smim(i-1, j-1) < 0) || (smim(i+1, j+1) < 0 && smim(i-1, j-1) >= 0)
                zc(i, j) = smim(i, j);
            elseif (smim(i-1, j+1) >= 0 && smim(i+1, j-1) < 0) || (smim(i-1, j+1) < 0 && smim(i+1, j-1) >= 0)
                zc(i, j) = smim(i, j);
            end
        end
    end
end
otpt=uint8(zc);
otptth= otpt>105;
axes(handles.ImgT);
imshow(otptth);

% --------------------------------------------------------------------
function MODELE_ELECTROSTATIQUE_Callback(hObject, eventdata, handles)
% hObject    handle to MODELE_ELECTROSTATIQUE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
if size(img, 3) == 3
    img = rgb2gray(img);
end
lambda = 0.04;
sigma = 1; 
seuil = 100; 
r = 6; 
w = 5 * sigma;
[m, n] = size(img); 
imd = double(img);
dya = [-sqrt(2)/4 0 sqrt(2)/4; -1 0 1; -sqrt(2)/4 0 sqrt(2)/4];
dxa = dya'; % derive verticale
g = fspecial('gaussian', max(1, fix(5 * sigma)), sigma); % gaussien
Ixa = conv2(imd, dxa, 'same');
Iya = conv2(imd, dya, 'same');
Ixa2 = conv2(Ixa.^2, g, 'same');
Iya2 = conv2(Iya.^2, g, 'same');
Ixya = conv2(Ixa.*Iya, g, 'same');
R = Ixa2 .* Iya2 - Ixya.^2 - lambda * (Ixa2 + Iya2).^2;

R1 = (1000 / (max(max(R)))) * R; %normalisation
[u, v] = find(R1 <= seuil);
nb = length(u);

for k = 1:nb
    R1(u(k), v(k)) = 0;
end

R11 = zeros(m + 2 * r, n + 2 * r);
R11(r + 1:m + r, r + 1:n + r) = R1;

[m1, n1] = size(R11);
for i = r + 1:m1 - r
    for j = r + 1:n1 - r
        fenetre = R11(i - r:i + r, j - r:j + r);
        ma = max(max(fenetre));
        if fenetre(r + 1, r + 1) < ma
            R11(i, j) = 0;
        end
    end
end

axes(handles.ImgT);
imshow(img);
hold on;
R11 = R11(r + 1:m + r, r + 1:n + r);
[x, y] = find(R11);
nb = length(x);
plot(y, x, '.r')


% --------------------------------------------------------------------
function Detection_ligne_Callback(hObject, eventdata, handles)
% hObject    handle to Detection_ligne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
imageIn=handles.courant_data;
if length( size(imageIn,3))>=3
    imageIn=rgb2gray(imageIn);
end
BW = edge(imageIn,'canny');
[H,theta,rho] = hough(BW);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,theta,rho,P,'FillGap',5,'MinLength',7);
axes(handles.ImgT);
imshow(imageIn), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red')


% --------------------------------------------------------------------
function Detection_cercle_Callback(hObject, eventdata, handles)
% hObject    handle to Detection_cercle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn = handles.courant_data;
if size(imageIn, 3) >= 3
    imageInGray = rgb2gray(imageIn);
else
    imageInGray = imageIn;
end
imageInGray = double(imageInGray) / 255;
e = edge(imageInGray, 'Canny');
% Trouvez les cercles dans l'image
%explication des parametres la recherche des cercles plus claire  que
%arriere plan avec un sensibilite de 0.9
[centers, radii] = imfindcircles(e, [15, 40], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);
axes(handles.ImgT);
imshow(imageIn), hold on
viscircles(centers, radii, 'EdgeColor', 'g');
hold off


% --------------------------------------------------------------------
function Homomorphique_Callback(hObject, eventdata, handles)
% hObject    handle to Homomorphique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn=handles.courant_data;
if size(imageIn, 3) == 3
        imageIn = rgb2gray(imageIn);
end
cim = double(imageIn);
cim = cim + 1;
lim = log(cim);
fim = fft2(lim);
lowg = 0.9; 
highg = 1.1; 
him = (highg - lowg) * (1 - exp(-0.002 * abs(fim))) + lowg;
him = him .* fim;
ifim = ifft2(him);
eim = exp(ifim);
axes(handles.ImgT);
imshow(uint8(eim))
