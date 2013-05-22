m = 30 # nombre des sous-variables

Contrainte01 = []
FonctionObject01 = []
Contrainte01.append('lambda x: True')
FonctionObject01.append('lambda x:x[0]')
f = lambda x:1+9*sum(map(lambda a:a/(m-1),x[1:]))
g = lambda x:(1-sqrt(x[0]/f(x) ))
FonctionObject01.append('lambda x:f(x)*g(x)')



Contrainte02 = []
FonctionObject02 = []
Contrainte02.append('lambda x: True')
FonctionObject02.append('lambda x:x[0]')
g2 = lambda x:(1-pow((x[0]/(f(x))),2))
FonctionObject02.append('lambda x:f(x)*g2(x)')

Contrainte03 = []
FonctionObject03 = []
Contrainte03.append('lambda x: True')
FonctionObject03.append('lambda x:x[0]')
g3 = lambda x:1-sqrt(x[0]/(f(x)))-(x[0]/(f(x)))*sin(10*pi*x[0])
FonctionObject03.append('lambda x:f(x)*g3(x)')



Contrainte04 = []
FonctionObject04 = []
Contrainte04.append('lambda x: True')
FonctionObject04.append('lambda x:x[0]')
FonctionObject04.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte05 = []
FonctionObject05 = []
Contrainte05.append('lambda x: True')
FonctionObject05.append('lambda x:x[0]')
FonctionObject05.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte06 = []
FonctionObject06 = []
Contrainte06.append('lambda x: True')
FonctionObject06.append('lambda x:1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6)')
FonctionObject06.append('lambda x:(1+9*pow(sum(x[1:])/(m-1),0.25))*(1-pow( ( 1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6) )/( 1+9*pow(sum(x[1:])/(m-1),0.25 )) ,2))')


