%%
%% This is file `WileySTAT-V1.cls'.
%% The original source files were from article.cls and Wiley Authoring template.
%%
%% Copyright (c) 2017 SPi Technologies.
%% For any queries in usage, mail us to LaTeXSupport@spi-global.com
%%
%% This file was generated from file(s) of the LaTeX base system.
%% --------------------------------------------------------------
%%
%% It may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License, either version 1.3c
%% of this license or (at your option) any later version.
%% The latest version of this license is in
%%    http://www.latex-project.org/lppl.txt
%% and version 1.3c or later is part of all distributions of LaTeX
%% version 2005/12/01 or later.

\def\update{v1.00}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{WileySTAT-V1}[2017/06/22 \update]

\newif\if@cmbrightfont
\DeclareOption{cmbright}{\@cmbrightfonttrue}

\newif\if@doublespace
\DeclareOption{doublespace}{\@doublespacetrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
\if@cmbrightfont
\usepackage[T1]{fontenc}
\RequirePackage[standard-baselineskips]{cmbright}
\fi

\if@doublespace
\RequirePackage[onehalfspacing]{setspace}
\fi

%% Load required packages
\RequirePackage{graphicx}
\RequirePackage[figuresright]{rotating}
\RequirePackage{pifont,latexsym,textcomp,ifthen,calc,textcase,booktabs,color}
\RequirePackage{amsfonts,amssymb,amsbsy,amsmath,amsthm}
\RequirePackage[errorshow]{tracefnt}
\RequirePackage{enumerate}%
\RequirePackage{algorithm, algorithmicx, algpseudocode}
\RequirePackage{listings}%
\RequirePackage{alltt}%
\RequirePackage{stfloats}%

%\RequirePackage{bm} %% For Bold math varients
\RequirePackage[errorshow]{tracefnt}

%% Followed as per Previous Template versions
\@twosidetrue
\flushbottom
\frenchspacing
\textwidth 175mm
\textheight 205mm
%\textheight 682pt
%\setlength\columnsep{24pt}

%Trim sizes
\setlength\voffset{-1in}
\setlength\hoffset{-1in}
\topmargin 4mm
\setlength\oddsidemargin{18mm}%back margin on odd pages
\setlength\evensidemargin{21.5mm}%fore margin on even pages
\setlength\paperwidth{210mm}
\setlength\paperheight{279mm}
%Needed to set PDF page size
%\special{papersize=210mm,279mm}

\parskip 6pt%\z@
\parindent 1em
\headheight 50pt
%\headsep    30pt
\headsep    24pt
\footskip 12mm

\brokenpenalty=10000%
\clubpenalty=9999%
\widowpenalty=9999%
\hyphenpenalty=9999%
\finalhyphendemerits=14000%
\doublehyphendemerits=75000%
\lefthyphenmin=3%
\righthyphenmin=3%
\pretolerance=1000%
\tolerance=2000%
\hfuzz=0.99pt%
\hbadness=8000
\vbadness=9000

\thinmuskip = 3mu
\medmuskip = 4mu
\thickmuskip = 5mu

\setcounter{topnumber}{10}
\def\topfraction{1}
\setcounter{bottomnumber}{10}
\def\bottomfraction{0.8}
\setcounter{totalnumber}{10}
\def\textfraction{0}
\renewcommand{\floatpagefraction}{0.95}
\setcounter{dbltopnumber}{10}
\renewcommand{\dblfloatpagefraction}{0.95}
\renewcommand{\dbltopfraction}{1}

\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xpt{13pt plus 0.10pt minus 0.10pt}%
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 6\p@ \@plus2\p@ \@minus2\p@
               \parsep 3\p@ \@plus2\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\normalsize

\renewcommand\small{%
   \@setfontsize\small\@viiipt{10}%
   \abovedisplayskip 8\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\renewcommand{\footnotesize}{\fontsize{9pt}{11pt}\selectfont}
%\renewcommand{\footnotesize}{\fontsize{7.5pt}{9.5pt}\selectfont}
\renewcommand\scriptsize{\@setfontsize\scriptsize\@viipt{9.5}}
\renewcommand\tiny{\@setfontsize\tiny\@vpt\@viipt}
\renewcommand\large{\@setfontsize\large\@xipt{13}}
\renewcommand\Large{\@setfontsize\Large\@xivpt{18}}
\renewcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\renewcommand\huge{\fontsize{22pt}{24pt}\selectfont}
\renewcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}

\newcommand{\sectionsize}{\fontsize{18pt}{20pt}\selectfont}
\newcommand{\subsectionsize}{\fontsize{15pt}{17pt}\selectfont}
\newcommand{\titlesize}{\fontsize{24pt}{26pt}\selectfont}
\newcommand{\authorsize}{\fontsize{15pt}{17pt}\selectfont}

\def\journalname{STAT}
\def\journalabb{stat}

\newbox\absbox
\def\abstract{\lrbox\absbox\minipage{\textwidth}%
  \normalsize{\bfseries Abstract:\hspace{0.75em}}%
  }
\def\endabstract{\copyrightline\endminipage\endlrbox}

\def\keywords#1{%
  \gdef\@keywords{{\sffamily%
    \textbf{Keywords:}}\hspace{0.75em}\parbox[t]{32pc}{#1}}}
\let\@keywords\@empty

\def\@received{??Please supply received date??}

\def\received#1{%
\gdef\@received{Received #1}}
%\let\@received\@empty
%\def\revised#1{%
%\gdef\@revised{Revised #1}}
%\let\@revised\@empty

\def\accepted#1{%
\gdef\@accepted{Accepted #1}}
\let\@accepted\@empty

\skip\footins 10pt plus  8pt
%\skip\footins 22pt plus  8pt
\gdef\footnoterule{}

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\long\def\@makefntext#1{\parindent 1em%
\noindent{$\m@th^{\@thefnmark}$}#1}

\def\corremail#1{%
  \gdef\@corraddr{%
    \footnotetext[1]{\bf Email: #1}\stepcounter{footnote}}}
\let\@corraddr\@empty
\def\corrauth{\footnotemark[1]}

%\def\corraddr#1{%
%  \gdef\@corraddr{%
%    \footnotetext[1]{\bf Email: #1}\stepcounter{footnote}}}
%\let\@corraddr\@empty
%\def\corrauth{\footnotemark[1]}

\def\email#1{%
  \gdef\@email{%
    \footnotetext[1]{{Email: #1}}\stepcounter{footnote}}}
\let\@email\@empty
\def\email{\footnotemark[2]}

\def\address#1{%
  \gdef\@address{%
    \footnotetext[0]{%\noindent{\large\hrule}%
    \rule{6pc}{.5pt}
    \vskip 4pt\noindent\bf #1}}}
\let\@address\@empty

%% Added articlenote
\def\articlenote#1{%
  \gdef\@articlenote{%
    \footnotetext[0]{\vskip -3pt\noindent #1}}}
\let\@articlenote\@empty

%\def\cgs#1{%
%  \gdef\@cgs{%
%    \footnotetext[0]{\par\noindent\itshape Contract/grant sponsor:  #1}}}
%\let\@cgs\@empty

\def\affilnum#1{${}^{\text{\bf#1}}$}
\def\affil#1{${}^{\text{\bf#1}}$}
\def\comma{$^{\text{\bf,}}$}

\renewcommand\maketitle{\par
  \begingroup
     \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
       \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{title}\label{FirstPage}\@address\@corraddr\@articlenote%\@cgs
  \endgroup
  \setcounter{footnote}{0}%
  \gdef\thefootnote{\arabic{footnote}}
  \global\let\address\relax
  \global\let\articlenote\relax
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\@address\@empty
  \global\let\@articlenote\@empty
  \global\let\corraddr\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{\vspace*{-30pt}%{\large\hrule}
%\vspace*{-20pt}%
%\noindent{\sffamily}\hfill\mbox{}\par
%\vspace{10pt}
%\noindent{\sffamily\normalsize(www.interscience.wiley.com) DOI: 10.1002/\DOI}\\
%\noindent{\sffamily\normalsize MOS subject classification: \@MOSsc}%
  \null
  \vskip 12pt %
    {\titlesize\raggedright\bfseries\@title \par}%
    \vskip 1.5em%
    {\authorsize
      \lineskip .5em%
      %
        \raggedright\sffamily\bfseries\@author
      \par}%
%            \vskip 12pt
%       \noindent\large \@received;\ \@accepted
 %{\normalsize
%      \lineskip .5em%
%      %
%        {\raggedright\emph\@address}
%      \par}%
    \vskip 14pt%
      {\noindent\normalsize\@keywords}
    {\lineskip 1.5em%
      %
    \vskip 2pt%
    {\noindent\usebox\absbox\par}
      \par}
      \vspace{1pt}
    \noindent%{\large\hrule}%\rule{\textwidth}{5.0pt}
    \par
  \vskip 2pc plus 12pt minus 6pt
  }

\def\startpage{\pageref{FirstPage}}
\def\endpage{\pageref{LastPage}}
\def\volumeyear{0000}
\def\volumenumber{00}
\def\DOI{stat.0000}

\gdef\copyrightline{Copyright \copyright\ \volumeyear\ John Wiley \& Sons, Ltd.}
\def\runninghead#1#2{\markboth{#1}{#2}}%{\color{blue}{Wiley StatsRef}}}

%\gdef\MOS#1{\gdef\@MOSsc{#1}}

\gdef\rhtext{{\Large\sffamily\bfseries \color{blue}{Wiley StatsRef}}}

%\gdef\ssboxnum{\colorbox{blue}{\sffamily\normalsize\textcolor{white}{\thepage}}}
\gdef\ssboxnum{{\sffamily\normalsize{\bfseries\thepage}}}

\def\ps@title{%
\def\@oddhead{\parbox{\textwidth}{{\Large\sffamily\bfseries {\color{blue}{Wiley StatsRef:}}\\[2pt]
 Statistics Reference Online %\\ Dissemination of Statistics Research
}\hfill{DOI: 10.1002/\DOI}\\[-2pt]
\rule{\textwidth}{2.0pt}
\noindent\normalfont%{\large\hrule}%\rule{\textwidth}{5.0pt}%
}}%
\let\@evenhead\@oddhead
\def\@oddfoot{\parbox[t]{\textwidth}{%
\noindent%{\large\hrule}%\rule{\textwidth}{5.0pt}
\par\vskip -2pt
{\sffamily\small{Wiley StatsRef}\ {\bfseries\volumeyear}, \volumenumber\ \startpage--\endpage\hfill\ssboxnum\hfill\copyrightline\\
\sffamily\textit{Prepared using \textsf{WileySTAT.cls} [Version: \update]}\mbox{}\\
}}}
\let\@evenfoot\@oddfoot
}

\def\ps@sim{%
\def\@evenhead{\parbox{\textwidth}{{\sffamily\rhtext\hfill\normalsize\leftmark}\\[-2pt]
\noindent{\large\hrule}%\rule{\textwidth}{5.0pt}
}}
\def\@oddhead{\parbox{\textwidth}{{\sffamily\normalsize\rightmark\hfill\rhtext}\\[-2pt]
\noindent{\large\hrule}%\rule{\textwidth}{5.0pt}
}}
\def\@evenfoot{\parbox[t]{\textwidth}{%\noindent{\large\hrule}%\rule{\textwidth}{5.0pt}
%\par
\sffamily\copyrightline\hfill\ssboxnum\hfill{Wiley StatsRef}\ {\bfseries\volumeyear}, \volumenumber\ \startpage--\endpage\\
\sffamily\textit{Prepared using \textsf{WileySTAT.cls}}}}
\def\@oddfoot{\parbox[t]{\textwidth}{%\noindent{\large\hrule}%\rule{\textwidth}{5.0pt}
%\par
\sffamily\small {Wiley StatsRef}\ {\bfseries\volumeyear}, \volumenumber\ \startpage--\endpage\hfill\ssboxnum\hfill\copyrightline\\
\sffamily\textit{Prepared using \textsf{WileySTAT.cls}}}}}

\let\@hangfrom\noindent %Remove indent on section turnovers
\renewcommand{\@seccntformat}[1]{{\csname the#1\endcsname.}\hspace{0.5em}}

\newdimen\@bls
\@bls=\baselineskip

\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {8pt plus 0pt minus 0pt}%
                                   {.3\@bls\@afterindentfalse}%
                                   {\normalfont\sectionsize\raggedright\sffamily\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {6pt plus 0pt minus 0pt}%
                                     {2pt\@afterindentfalse}%
                                     {\normalfont\subsectionsize\raggedright\sffamily\itshape}}

%\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
%                                     {\@bls plus .3\@bls minus .1\@bls}%
%                                     {6pt\@afterindentfalse}%
%                                     {\normalfont\normalsize\raggedright\sffamily\itshape}}

\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {\@bls plus .2\@bls}%
                                     {1pt\@afterindentfalse}%
                                     {\normalfont\normalsize\sffamily\itshape}}

\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                     {.5\@bls plus .2\@bls}%
                                     {-1em}%
                                     {\normalfont\normalsize\bfseries}}

\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                      {.25\@bls plus .2\@bls}%
                                      {-1em}%
                                      {\normalfont\normalsize\bfseries}}

\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\itemsep 0pt\parsep 0pt
         \def\makelabel##1{\hss\llap{##1}}}\fi}

\let\endenumerate =\endlist

\def\itemize{\ifnum \@itemdepth >3 \@toodeep\else \advance\@itemdepth \@ne
\edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
\list{\csname\@itemitem\endcsname}{\itemsep 0pt\parsep 0pt
\def\makelabel##1{\hss\llap{##1}}}\fi}

\let\enditemize =\endlist

\usepackage[tableposition=top]{caption}
\DeclareCaptionLabelSeparator{jwperiod}{.\hspace*{0.5ex}}
\captionsetup[figure]{font=footnotesize,labelfont=bf,labelsep=jwperiod,justification=justified,singlelinecheck=true}
\captionsetup[table]{font=footnotesize,labelfont=bf,labelsep=jwperiod,justification=centerlast,singlelinecheck=true}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theorem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\RequirePackage{amsthm}
\def\thm@space@setup{%
  \thm@preskip=7\p@ \thm@postskip=7\p@
}
\thm@headpunct{.}

\newtheoremstyle{plain}% Numbered
{7pt}% Space above
{7pt}% Space below
{\normalfont\itshape}% Body font
{\parindent}% Indent amount
{\bfseries}% Theorem head font
{.\quad}% Punctuation after theorem head
{0pt}% Space after theorem headi
{}% Theorem head spec (can be left empty, meaning `normal')

\renewenvironment{proof}[1][\proofname]{%\par
  \pushQED{\qed}%
  \normalfont \topsep3\p@\@plus2\p@\relax
  {\bfseries #1}.\quad\itshape
}{%
  \popQED\endtrivlist\@endpefalse
}

%%%%%%%%%%%%%

\newcommand\ack{\section*{Acknowledgement}}
\newcommand\acks{\section*{Acknowledgements}}

\newcommand\BibTeX{{\rmfamily B\kern-.05em \textsc{i\kern-.025em b}\kern-.08em T\kern-.1667em\lower.7ex\hbox{E}\kern-.125emX}}
\renewcommand\refname{References}

\renewenvironment{thebibliography}[1]{%
\section*{\refname}
\list{{\arabic{enumi}}}{\def\makelabel##1{\hss{##1.}}\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=0\p@
\labelsep=1ex\itemindent=0\p@
\settowidth\labelwidth{\small[#1]}%
\leftmargin\labelwidth
\advance\leftmargin\labelsep
\advance\leftmargin -\itemindent
\usecounter{enumi}}\small
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}{\endlist}

\def\biog{\section*{Author's Biography}\small}
\def\biogs{\section*{Authors' Biographies}\small}

\AtEndDocument{\label{LastPage}}

\hyphenation{com-mu-ni-ca-tions}

\pagestyle{sim}
\normalsize
\sloppy

%% Appendix
\renewcommand\appendix{%\par
    \setcounter{section}{0}%
    \setcounter{equation}{0}%
    \setcounter{table}{0}%
    \setcounter{figure}{0}%
    \gdef\thesection{\@Alph\c@section}
    \@addtoreset{equation}{section}
    \gdef\theequation{\@Alph\c@section\arabic{equation}}
    \@addtoreset{table}{section}
    \gdef\thetable{\@Alph\c@section\arabic{table}}
    \@addtoreset{figure}{section}
    \gdef\thefigure{\@Alph\c@section\arabic{figure}}
  }

%% Figure source
\newcommand{\figsource}[2][\textit{Source:}]{\null\strut\par#1\enskip#2}
\arraycolsep2pt

\endinput
%% Version 1.00 -- Initial updates
