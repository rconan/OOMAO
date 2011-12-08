\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{RSAA SDN}[2010/12/15 1.1]


\RequirePackage{ifthen}
\RequirePackage[dvips]{graphicx}
\RequirePackage{color}
\usepackage{amsmath,amssymb}
\usepackage{lastpage}
\usepackage{booktabs}
\usepackage{tikz}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[twoside,a4paper,10pt]{article}

\newcommand{\Project}[1]  {\renewcommand{\DH@project}{#1}}
\newcommand{\DocNumber}[1]  {\renewcommand{\DH@number}{#1}}
\newcommand{\DocVersion}[1] {\renewcommand{\DH@version}{#1}}
\newcommand{\DH@number}{} 
\newcommand{\DH@version}{}
\newcommand{\DH@project}{}
\@onlypreamble{\DocNumber}
\@onlypreamble{\DocVersion}
\@onlypreamble{\Project}

\newcommand{\anu}{\textsc{The Australian National University}}
\newcommand{\logo}{
\ifthenelse{\value{page}=1}%
{
  \begin{minipage}{0.5\linewidth}
    \begin{tikzpicture}
\draw[very thick] (0,0) circle (0.2cm);
\draw[very thick] (0,0) circle (0.75mm);
\foreach \x in {0,60,120,180,240,300}
\draw[very thick] (\x:0.45cm) circle (0.2cm);
    \node[anchor=west] at (0.65,0) {\textbf{\Huge{\DH@project}}};
    \end{tikzpicture}
  \end{minipage}
}%
{\textbf{\huge{\DH@project}}}
}

\newcommand{\docinfo}{%
  \begin{tabular}{r}
%    System Design Note\\
    \DH@number \\
    Version \DH@version
  \end{tabular}
}
\newcommand{\anudocinfo}{%
%   \begin{tabular}{r@{}}
%     \anu\\\vspace{1ex}
% %    System Design Note\\
% %    \DH@number \\
%     Version \DH@version
%   \end{tabular}
\includegraphics[scale=0.5]{anu_logo.jpg}
}

\newcommand{\oddinfo}{
\ifthenelse{\value{page}=1}%
{
    \anudocinfo
}{%
  \begin{tabular}{r@{}}
%    System Design Note\\
    \DH@number \\
    Version \DH@version
  \end{tabular}
}}

\newcommand{\eveninfo}{
  \begin{tabular}{@{}l}
%    System Design Note\\
    \DH@number \\
    Version \DH@version
  \end{tabular}
}

%\addtolength{\hoffset}{-0.5cm} 
\addtolength{\textwidth}{3cm}
\addtolength{\headheight}{5mm}
\addtolength{\evensidemargin}{-2cm}
\addtolength{\oddsidemargin}{-1cm}

 \renewcommand\maketitle{%
% 
   \begin{minipage}{1.0\linewidth}
     \vspace{20ex}
     \begin{center}
               {\setlength{\baselineskip}%
{1.5\baselineskip}
\textbf{\Huge {\@title}}\par}
             \end{center}\vspace{10ex}
             \begin{center}
               {\large Research School of Astronomy and Astrophysics\\
               ANU College of Physical and Mathematical Sciences\\
               The Australian National University\\}
             \end{center}
     \vspace{5ex}
     \begin{center}
       \DH@number
     \end{center}
     \vspace{8ex}
\end{minipage} }


 \newcommand{\ps@dochead}{%
 \renewcommand{\@evenhead}{%
   \begin{minipage}{\textwidth}
     \begin{tabular*}{\textwidth}{@{}l@{\extracolsep{\fill}}r@{}}%
       \eveninfo     & \logo \\
      \end{tabular*}\vspace{2ex} \rule{\textwidth}{1pt}
   \end{minipage}}
 \renewcommand{\@oddhead}{%
   \begin{minipage}{\textwidth}
     \begin{tabular*}{\textwidth}{@{}r@{\extracolsep{\fill}}l@{}}%
       \logo    &    \oddinfo     \\
      \end{tabular*}\vspace{2ex} \ifthenelse{\value{page}=1}{\rule{\textwidth}{2mm}\vspace{4ex}}{\rule{\textwidth}{1pt}} %
   \end{minipage}}
 \renewcommand{\@oddfoot}{
   \begin{minipage}{1.0\linewidth}
     \begin{center}
       Page \thepage\ of \pageref{LastPage}
     \end{center}
   \end{minipage}
}
 \renewcommand{\@evenfoot}{
   \begin{minipage}{1.0\linewidth}
     \begin{center}
       Page \thepage\ of \pageref{LastPage}
     \end{center}
   \end{minipage}
}
 }
 \pagestyle{dochead}

 %\setlength{\headheight}{\baselineskip}
 %\setlength{\headsep}{3em}

%  \def\tableofcontents{
%    \section*{Contents}
%    \@starttoc{toc}
% %   \clearpage
%  }

 \newenvironment{revision}
{\vspace{2ex}\section*{\centering Revision History}%
\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}|l|p{0.2\textwidth}|l|p{0.3\textwidth}|}%
\hline
\textbf{Version No.} & \textbf{Author \& Date} & \textbf{Approval \& Date} & \textbf{Description} \\\hline}
{\hline
\end{tabular*}}

 \newenvironment{documents}
{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}|l|l|p{0.3\textwidth}|l|l|}%
\hline
\textbf{Document ID} & \textbf{Source} & \textbf{Title} & \textbf{Version} & \textbf{Date} \\\hline}
{\hline
\end{tabular*}}
