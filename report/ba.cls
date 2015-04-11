%%
%% This is file `ba.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% ba.dtx  (with options: `bacls')
%% ===========================================================================
%% Bayesian Analysis Journal document class and packages for use with LaTeX 2e
%%
%% STYLE FILES CREATED FOR AND BY THE STAT JOURNAL 
%% 
%% Some macros and environments are modified versions of those given in
%% classes.dtx v1.4a 1999/01/07, copyrighted by The LaTeX3 Project.
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
\NeedsTeXFormat{LaTeX2e}[1996/06/01]
\ProvidesClass{ba}
              [2001/12/17 v1.09
 Bayesian Analysis Journal document class]
\def\ba@Warning#1{\ClassWarning{ba}{#1}}
\def\ba@WarningNoLine#1{\ClassWarningNoLine{ba}{#1}}
\def\ba@Info#1{\typeout{#1}\ClassInfo{ba}{#1}}
\newif\ifba@opt@debug
\DeclareOption{debug}{\ba@opt@debugtrue}
\DeclareOption{nodebug}{\ba@opt@debugtrue}
\newif\ifba@opt@book
\DeclareOption{book}{\ba@Info{book is base class}\ba@opt@booktrue}
\DeclareOption{nobook}{\ba@Info{article is base class}\ba@opt@bookfalse}
\newif\ifba@opt@bib
\ba@opt@bibtrue
\DeclareOption{bib}{\relax}
\DeclareOption{nobib}{\ba@opt@bibfalse}
\DeclareRobustCommand\em
        {\@nomath\em \ifdim \fontdimen\@ne\font >\z@
                       \upshape \else \slshape \fi}
\newif\ifba@opt@fancy
\DeclareOption{fancy}{\ba@opt@fancytrue}
\DeclareOption{nofancy}{\ba@opt@fancyfalse}
\newif\ifba@opt@crop
\DeclareOption{crop}{\ba@opt@croptrue}
\DeclareOption{nocrop}{\ba@opt@cropfalse}
\newif\ifba@opt@reprints
\DeclareOption{reprints}{\ba@opt@reprintstrue}
\DeclareOption{noreprints}{\ba@opt@reprintsfalse}
\newif\ifba@pdftex
\ifx\pdfoutput\ba@undefined
\else
  \ba@pdftextrue
\fi
\def\ba@opt@ptsize{10pt}
\DeclareOption{10pt}{\relax}
\DeclareOption{11pt}{\def\ba@opt@ptsize{11pt}}
\DeclareOption{12pt}{\def\ba@opt@ptsize{12pt}}
\let\baefinalopt=\ba@undefined
\let\bafinalopt=\ba@undefined
\DeclareOption{default}{
  \ba@opt@booktrue
  \ba@opt@croptrue
  \let\baefinalopt=\ba@undefined
  \let\bafinalopt=\ba@undefined
}
\DeclareOption{final}{
  \ba@opt@croptrue
  \ba@opt@fancytrue
  \let\baefinalopt=\ba@undefined
  \let\bafinalopt=\relax
}
\DeclareOption{efinal}{
  \ba@opt@croptrue
  \ba@opt@fancytrue
  \let\baefinalopt=\relax
  \let\bafinalopt=\ba@undefined
}
\DeclareOption*{\ba@Warning{{\CurrentOption} is not a valid option}}
\ExecuteOptions{default}
\ProcessOptions\relax
\ifba@opt@book
  \PassOptionsToClass{\ba@opt@ptsize}{book}
  \LoadClass{book}
\else
  \PassOptionsToClass{\ba@opt@ptsize,twoside}{article}
  \LoadClass{article}
\fi
\RequirePackage{calc}[1998/07/07]
\RequirePackage{ifthen}[1999/09/10]
\ifba@opt@fancy
  \RequirePackage{fancyhdr}% no date
\fi
\ifba@pdftex
  \ba@Warning{crop.sty not being used with pdflatex}
  \ba@opt@cropfalse
\fi
\ifba@opt@crop
  \RequirePackage[lettercenter]{crop}[1999/05/15]
  \renewcommand*\CROP@@info{\relax}
\fi
\ifba@opt@bib
  \RequirePackage{natbib}[1999/05/28]
  \RequirePackage{chapterbib}[1999/04/29]
  \bibpunct{(}{)}{;}{a}{}{,}
\fi
\ifba@opt@book
\else
  \newcounter{chapter}
\fi
\renewcommand{\thepart}{\arabic{part}}
\renewcommand{\thechapter}{\arabic{chapter}}
\renewcommand{\thesection}{\arabic{section}}
\renewcommand{\thesubsection}{\thesection.\arabic{subsection}}
\renewcommand{\theequation}{\arabic{equation}}
\renewcommand{\thefigure}{\arabic{figure}}
\renewcommand{\thetable}{\arabic{table}}
\pagestyle{headings}
\def\ps@headings{%
  \ifba@opt@book
  \else
    \newcommand{\chaptermark}[1]{\markboth{\textbf{##1}}{}}
  \fi
  \renewcommand{\chaptermark}[1]{\markboth{\textbf{\protect\short@title}}{}}
  \renewcommand{\sectionmark}[1]{\markright{\protect\ba@short@author}}
  \renewcommand{\subsectionmark}[1]{\markright{\protect\ba@short@author}}
  \renewcommand{\subsubsectionmark}[1]{\markright{\protect\ba@short@author}}
}
\newif\ifba@anyinserts
\def\inserttag{notag}
\def\ba@tagfilename{\jobname.tag}
\newwrite\ba@tagfile
\immediate\openout\ba@tagfile\ba@tagfilename
\def\ba@inserttypw@previous{\relax}
\newcommand{\inserttype}[2][notag]{%
  \clearpage
  \ifthenelse{\equal{#1}{notag}}{%
    \def\inserttag{notag\thechapter}%
  }{
    \def\inserttag{#1}%
  }
  \let\ba@inserttype@previous\ba@inserttype
    \def\ba@inserttype{article}%
    \let\@maketitle\ba@maketitle@art
    \babibsection
  \ifthenelse{\equal{#2}{column}}{%
    \def\ba@inserttype{column}%
    \let\@maketitle\ba@maketitle@col
    \babibsection
  }{}
  \ifthenelse{\equal{#2}{note}}{%
    \def\ba@inserttype{note}%
    \let\@maketitle\ba@maketitle@not
    \babibsection[star]
  }{}
  \ifthenelse{\equal{#2}{comment}}{%
    \def\ba@inserttype{comment}%
    \let\@maketitle\ba@maketitle@com
    \babibsection[star]
  }{}%
  \ifthenelse{\equal{#2}{update}}{%
    \def\ba@inserttype{update}%
    \let\@maketitle\ba@maketitle@upd
    \babibsection[star]
    \author{\relax}{\relax} % updates do not have specific authors
  }{}%
}
\newcommand{\ba@maketitle@art}{%
  \let\and\cand
  {\parindent 0pt \centering \normalfont\sffamily
   \LARGE\bfseries\ba@placetitle\par\nobreak
   \vskip 16pt
  }
  {\parindent 0pt \centering \normalfont%
   \begin{tabular}[t]{c}
   \@author
   \end{tabular}
   \par\nobreak%
  }%
}
\newcommand{\columnfig}[1]{\def\@columnfig{#1}}
\let\emptycolumnfig\relax
\let\@columnfig\emptycolumnfig
\newbox\columnpicbox
\newcount\columnpichght
\newcommand{\ba@maketitle@col}{%
  \let\and\cand
  \ifx\@columnfig\emptycolumnfig
    \ba@maketitle@art
  \else
    \setbox\columnpicbox=\hbox{\epsfig{file=\@columnfig}}
    \columnpichght=\ht\columnpicbox
    {\parindent 0pt \centering \normalfont\sffamily
     \LARGE\bfseries \ba@placetitle\par\nobreak
     \vskip 16pt
    }
    {\parindent 0pt
     \normalfont
     \begin{minipage}[b][\the\columnpichght sp][c]{\textwidth/2-1em}
     \raggedleft\@author
     \end{minipage}
     \begin{minipage}[b][\the\columnpichght sp][c]{\textwidth/2-1em}
     \unhbox\columnpicbox
     \end{minipage}
     \par\nobreak
    }
    % reset column figure macro to be empty
    \let\@columnfig\emptycolumnfig
  \fi
  \vskip -6pt
}
\newcommand{\ba@maketitle@not}{%
  \let\and\rand
  \noindent
  {\raggedleft\normalfont\sffamily
   \Large\bfseries \ba@placetitle\par\nobreak
   \vskip 10pt
  }
  {\parindent 0pt \raggedleft\normalfont
   \begin{tabular}[t]{@{}r@{}}
   \@author
   \end{tabular}
   \par\nobreak
   \vskip -6pt
  }
}
\let\ba@maketitle@com=\ba@maketitle@not
\newcommand{\ba@maketitle@upd}{%
  \noindent
  {\raggedright\normalfont\sffamily
   \Large\bfseries \ba@placetitle\par\nobreak
   \vskip 20pt
  }
}
\newcommand{\ba@maketitle@def}{%
  {\parindent 0pt \centering \normalfont\sffamily
   \LARGE\bfseries \ba@placetitle\par\nobreak
   \vskip 16pt
  }
  {\parindent 0pt \centering \normalfont
   \begin{tabular}[t]{c}
   \@author
   \end{tabular}
   \par\nobreak
   \vskip -6pt
  }
}
\let\@maketitle\ba@maketitle@def
\let\ba@short@author@previous\relax
\def\ba@short@author{\relax}
\renewcommand{\author}[2]{%
  \let\ba@short@author@previous\ba@short@author
  \def\baauthor{#1}%
  \def\ba@short@author{#1}%
  \def\@author{#2}%
}
\def\land{%
  \end{tabular}%
  \hskip 1em \@plus.17fil%
  \begin{tabular}[t]{@{}l@{}}%
}
\def\cand{%
  \end{tabular}%
  \hskip 1em \@plus.17fil%
  \begin{tabular}[t]{@{}c@{}}%
}
\def\rand{%
  \end{tabular}%
  \hskip 1em \@plus.17fil%
  \begin{tabular}[t]{@{}r@{}}%
}
\renewcommand{\title}[2][empty]{%
  \gdef\batitle{#2}%
  \gdef\@title{#2}%
  \ifthenelse{\equal{#1}{empty}}{
    \gdef\short@title{#2}
  }{
    \gdef\short@title{#1}
  }%
}
\def\ba@placetitle{%
   \ifx\hypertarget\ba@undefined%
    \@title%
   \else%
    \hypertarget{\inserttag}{\@title}%
   \fi%
}
\renewcommand{\part}[2][]{%
  \ifthenelse{\equal{#1}{noclear}}{\relax}{\clearpage}%
  \refstepcounter{part}
  \ifnum \thepart>1
    \addtocontents{toc}{\protect\addvspace{1.5\baselineskip}}%
  \fi
  \addcontentsline{toc}{part}{\bfseries #2}
  \addtocontents{toc}{\protect\addvspace{\baselineskip}}%
}
\def\baaddcontentsline{%
  \addcontentsline{toc}{chapter}{%
    \@title\dotfill\newline%
    \hbox{\hskip20pt}\dotfill\ba@short@author%
  }%
}
\renewcommand{\@pnumwidth}{\@tocrmarg}
\newlength{\ba@len@max@title}
\setlength{\ba@len@max@title}{\textwidth}
\addtolength{\ba@len@max@title}{\@tocrmarg}
\addtolength{\ba@len@max@title}{-1em}
\newlength{\ba@len@a}
\newlength{\ba@len@b}
\renewcommand{\maketitle}{%
  \ba@pp@set
  \ba@post@reprints{\ba@short@author@previous}
  \clearpage\ba@pp@setfirst
  \bathisheadings%
    {\small Bayesian Analysis (\ba@year)\newline%
     \textbf{\ba@volume}, Number \ba@number \ba@pp@get}
    {}
    {}
    {\small\copyright{} \ba@year{} International Society for Bayesian Analysis}
    {}
    {\inserttag}
  \pagestyle{headings}%
  \refstepcounter{chapter}%
  \setcounter{section}{0}
  \setcounter{equation}{0}
  \setcounter{figure}{0}
  \setcounter{table}{0}
  \vskip\baselineskip\par%
  \@maketitle%
  \chaptermark{\short@title}%
  \sectionmark{\ba@short@author}%
  \immediate\write\ba@tagfile{insert \thechapter --\inserttag: \@title}
  \ba@Info{New \ba@inserttype: (chapter \thechapter: dest \inserttag)}%
  \ifthenelse{\equal{\ba@inserttype}{update}}{%
    \part[noclear]{\short@title}%
  }{%
    \settowidth{\ba@len@a}{\normalfont\@title}%
    \settowidth{\ba@len@b}{\normalfont\ba@short@author}%
    \addtolength{\ba@len@a}{\ba@len@b}%
    \ifdim\ba@len@a<\ba@len@max@title%
      \addcontentsline{toc}{chapter}{\@title\dotfill\ba@short@author}%
    \else
      \addtolength{\ba@len@a}{-\ba@len@b}%
      \ifdim\ba@len@a<\ba@len@max@title%
        \addcontentsline{toc}{chapter}{%
          \@title\dotfill\newline%
          \hbox{\hskip0pt}\dotfill\ba@short@author%
        }%
      \else
        \addtolength{\ba@len@a}{-\ba@len@max@title}
        \addtolength{\ba@len@a}{\ba@len@b + \parindent}
        \ifdim\ba@len@a<\ba@len@max@title%
          \addcontentsline{toc}{chapter}{\@title\dotfill\hbox{\ba@short@author}}%
        \else
          \addtolength{\ba@len@a}{-\ba@len@max@title}
          \ifdim\ba@len@a<\ba@len@max@title%
            \addcontentsline{toc}{chapter}{%
              \@title\dotfill\newline%
              \hbox{\hskip0pt}\dotfill\ba@short@author%
            }%
          \else
            \baaddcontentsline%
          \fi
        \fi
      \fi
    \fi
  }
  \ba@anyinsertstrue%
}
\ifba@opt@reprints
  \def\ba@reprints@filename{\jobname.reprints}
  \newwrite\ba@reprints@file
  \openout\ba@reprints@file\ba@reprints@filename
\fi
\def\ba@pp@basefile{\jobname.pp}
\newwrite\ba@pp@file
\newcommand{\ba@pp@get}{%
  \IfFileExists{\ba@pp@basefile\thechapter}{%
    \input{\ba@pp@basefile\thechapter}%
  }{}%
}
\newcounter{ba@pp@thispage}
\newcommand{\ba@post@reprints}[1]{%
  \ifba@opt@reprints
    \ifba@anyinserts
      \ifthenelse{%
        \equal{\ba@inserttype@previous}{article}%
      }{%
        \immediate\write\ba@reprints@file{}
        \immediate\write\ba@reprints@file{#1}
        \immediate\write\ba@reprints@file%
        {pp: \ba@pp@first--\theba@pp@thispage, copies 25}
      }{}
      \ifthenelse{%
        \equal{\ba@inserttype@previous}{column}%
      }{%
        \immediate\write\ba@reprints@file{}
        \immediate\write\ba@reprints@file{#1}
        \immediate\write\ba@reprints@file%
        {pp: \ba@pp@first--\theba@pp@thispage, copies 25}
      }{}
    \fi
  \fi
}
\newcommand{\ba@pp@set}{%
  \ifba@anyinserts
    \setcounter{ba@pp@thispage}{\thepage}
    \ifnum \theba@pp@thispage>1 \addtocounter{ba@pp@thispage}{-1} \fi
    \immediate\openout\ba@pp@file\ba@pp@basefile\thechapter
    \ifnum \theba@pp@thispage>\ba@pp@first
      \immediate\write\ba@pp@file{, pp. \ba@pp@first--\theba@pp@thispage}
    \else
      \immediate\write\ba@pp@file{, p. \theba@pp@thispage}
    \fi
    \immediate\closeout\ba@pp@file
  \fi
}
\newcommand{\ba@pp@setfirst}{%
  \edef\ba@pp@first{\thepage}
}
\def\ba@pp@first{0}
\AtEndDocument{
  \ba@pp@set
  \inserttype{update}
  \ba@post@reprints{\ba@short@author@previous}%
  \ifba@opt@reprints
    \closeout\ba@reprints@file
  \fi
}
\renewcommand{\section}{\@startsection%
  {section}%
  {1}%
  {0pt}%
  {-\baselineskip}%
  {0.25\baselineskip}%
  {\normalfont\Large\sffamily\bfseries}%
}
\let\ba@section=\section
\renewcommand{\subsection}{\@startsection%
  {subsection}%
  {2}%
  {0pt}%
  {-\baselineskip}%
  {0.25\baselineskip}%
  {\normalfont\large\sffamily\bfseries}%
}
\let\ba@subsection=\subsection
\renewcommand{\subsubsection}{\@startsection%
  {subsubsection}%
  {3}%
  {0pt}%
  {-\baselineskip}%
  {0.25\baselineskip}%
  {\normalfont\normalsize\sffamily\bfseries}%
}
\let\ba@subsubsection=\subsubsection
\renewcommand{\paragraph}{\@startsection%
  {paragraph}%
  {4}%
  {0pt}%
  {-\baselineskip}%
  {-0.25\baselineskip}%
  {\normalfont\normalsize\sffamily\bfseries}%
}
\let\ba@paragraph=\paragraph
\def\contentsname{Contents}
\renewcommand{\tableofcontents}[1][]{%
  \clearpage
  \bathisheadings%
    {\normalfont\Large Volume \ba@volume}
    {\normalfont\Large Number \ba@number}
    {\normalfont\Large\ba@year}
    {}
    {}
    {}
    {}
  \begingroup
    \parskip 0pt
    {\parindent 0pt \centering \normalfont
     \Huge\scshape Bayesian Analysis
     \vskip \baselineskip
    }
  {\parskip 2pt \@starttoc{toc}}%
  \endgroup
  \clearpage
  \thispagestyle{empty}
  #1
  \pagenumbering{arabic}
}
\setcounter{tocdepth}{0}
\renewcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{1em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       {\large \bfseries #1}\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
    \endgroup
  \fi
}
\ifba@opt@book\else
  \newcommand{\l@chapter}{\relax}
\fi
\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}
\ifba@opt@fancy
  \renewcommand{\headrulewidth}{0pt}
\fi
\newcommand{\bathisheadings}[6]{%
  \ifba@opt@fancy
    \thispagestyle{fancy}
    \ba@fancyheadings{#1}{#2}{#3}{#4}{#5}{#6}
  \else
    \thispagestyle{myheadings}
    \ba@myheadings{#1}{#2}{#3}{#4}{#5}{#6}
  \fi
}
\newcommand{\ba@headings}[6]{%
  \ifba@opt@fancy
    \ba@fancyheadings{#1}{#2}{#3}{#4}{#5}{#6}
  \else
    \ba@myheadings{#1}{#2}{#3}{#4}{#5}{#6}
  \fi
}
\newcommand{\ba@fancyheadings}[6]{%
  \lhead{#1}
  \chead{#2}
  \rhead{#3}
  \lfoot{#4}
  \cfoot{#5}
  \rfoot{#6}
}
\newcommand{\ba@myheadings}[6]{%
  \renewcommand{\ps@myheadings}{%
    \renewcommand{\@oddhead}{{#1\hfil#2\hfil#3}}
    \renewcommand{\@evenhead}{{#3\hfil#2\hfil#2}}
    \renewcommand{\@oddfoot}{{#4\hfil#5\hfil#6}}%
    \renewcommand{\@evenfoot}{{#6\hfil#5\hfil#4}}%
  }
}
\newcommand{\keywords}[1]{%
  \par%
  \noindent%
  {\small {\sffamily\bfseries Keywords:} #1}%
}
\newenvironment{summary}{%
  \small%
  \noindent%
  {\sffamily\bfseries Summary.}%
}{\par}
\ifba@opt@book
  \newenvironment{abstract}{\relax}{\relax}
\fi
\def\abstractname{Abstract}%
\renewenvironment{abstract}{%
  \small%
  \quotation%
  \parskip 0.5\baselineskip%
  {\noindent\sffamily\bfseries\abstractname.}
}
{ \endquotation%
  \par%
}
\newenvironment{acknowledgement}{%
  \vskip 0.5\baselineskip\small%
  {\noindent\normalfont\sffamily\bfseries Acknowledgments}\par%
  \begingroup\parindent 0pt\parskip 0.5\baselineskip%
}
{\endgroup}
\newenvironment{aboutauthor}{%
  \vskip 0.5\baselineskip\small%
  {\noindent \normalfont\sffamily\bfseries About the Author}\par%
  \begingroup\parindent 0pt\parskip 0.5\baselineskip%
}
{\endgroup}
\newenvironment{aboutauthors}{%
  \vskip 0.5\baselineskip\small%
  {\noindent\normalfont\sffamily\bfseries About the Authors}\par%
  \begingroup\parindent 0pt\parskip 0.5\baselineskip%
}
{\endgroup}
\AtBeginDocument{
  \let\ba@caption\caption
  \def\caption{%
    \ifthenelse{\equal{\@captype}{figure}}{\unskip}{\relax}%
    \ba@caption
  }
}
\ifba@opt@bib
  \newdimen\bibindent
  \setlength\bibindent{1.5em}
  \def\bibname{References}
  \newcommand\babibsection[1][]{%
    \ifthenelse{\equal{#1}{}}{%
      \def\bibsection{\section*{\bibname}}%
    }{%
      \def\bibsection{\section*{\bibname}}%
    }
  }
\fi
\def\defaultoffset{1in}
\newlength{\ba@pgheight}
\newlength{\ba@pgwidth}
\setlength{\ba@pgheight}{9.25in}
\setlength{\ba@pgwidth}{7.25in}
\def\ba@setpaperdim{
  \setlength{\paperheight}{\ba@pgheight}
  \setlength{\paperwidth}{\ba@pgwidth}
}
\newcommand{\bamakecrop}[1][baempty]{%
  \ifthenelse{\equal{#1}{baempty}}{%
    \ifx\baefinalopt\ba@undefined \crop[cam]%
    \else \crop[frame]%
    \fi%
  }{%
    \crop[#1]%
  }%
}
\ifba@opt@crop
  \ba@setpaperdim
\else
  \newcommand{\crop}[1][empty]{\ba@Info{crop.sty not loaded}}
\fi
\newlength{\ba@voffset}\setlength{\ba@voffset}{0.25in}
\newlength{\ba@hoffset}\setlength{\ba@hoffset}{\defaultoffset}
\ifba@pdftex
  \setlength{\paperheight}{11in}
  \setlength{\paperwidth}{8.5in}
  \setlength{\voffset}{0.5\paperheight-0.5\ba@pgheight}
  \setlength{\hoffset}{0.5\paperwidth-0.5\ba@pgwidth}
\fi
\setlength{\oddsidemargin}{0pt}
\setlength{\evensidemargin}{0pt}
\newlength{\ba@topmargin}\setlength{\ba@topmargin}{0pt}
\ifba@pdftex
  \setlength{\topmargin}{\ba@topmargin-\defaultoffset}
\else
  \setlength{\topmargin}{\ba@topmargin-\defaultoffset+\ba@voffset}
\fi
\setlength{\headheight}{44pt}
\setlength{\headsep}{14pt}
\setlength{\marginparsep}{0pt}
\setlength{\marginparwidth}{0pt}
\setlength{\footskip}{30pt}
\setlength{\textheight}%
  {\ba@pgheight-2\ba@voffset-\ba@topmargin-\headheight-\headsep-\footskip}
\setlength{\textwidth}%
 {\ba@pgwidth-2\ba@hoffset-\oddsidemargin-\marginparsep-\marginparwidth}
\ifba@opt@crop
  \relax
\else
  \ifba@pdftex
    \relax
  \else
    \setlength{\paperheight}{11in}
    \setlength{\paperwidth}{8.5in}
    \setlength{\topmargin}{0pt}
    \setlength{\voffset}{0pt}
    \setlength{\hoffset}{0.5\paperwidth-0.5\textwidth-1in}
  \fi
\fi
\newlength{\banoheadlength}
\setlength{\banoheadlength}{-\headheight}
\setlength{\parskip}{0.5\baselineskip}
\ifx\baefinalopt\ba@undefined \else
  \def\forcebftt{\relax}
\fi
\def\ba{Bayesian Analysis}
\def\tba{The {\ba}}
\newcommand{\basetissue}[4]{%
  \def\ba@volume{#1}
  \def\ba@number{#2}
  \def\ba@month{#3}
  \def\ba@year{#4}
}
\basetissue{1}{1}{October}{2004}
\AtBeginDocument{%
  \IfFileExists{editor.sty}{%
    \ba@Info{automatically inputting editor.sty}\usepackage{editor}%
  }{}%
  \IfFileExists{ba.cfg}{%
    \ba@Info{automatically inputting ba.cfg}\input{ba.cfg}%
  }{}%
}
% defs for Metadata
\def\fnms#1{#1}
\def\snm#1{#1}
\def\ead#1{}
\def\kwd#1{#1}
% 

\endinput
%%
%% End of file `ba.cls'.
