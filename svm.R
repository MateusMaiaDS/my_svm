#Setting training and and test toy example
set.seed(42)
split_ratio<-0.7
data(spam)
data<-mlbench::mlbench.circle(500,2) %>% as.data.frame
train_index<-sample(1:nrow(data),round(split_ratio*nrow(data)))
train<-data[train_index,]
test<-data[-train_index,]

library(kernlab)
mod_ksvm<-kernlab::ksvm(classes~.,data=train,C=1,kernel='vanilladot')
pred_ksvm<-predict(mod_ksvm,newdata=test) %>% as.vector()


#Calculating w
alfa_zero<-rep(0,nrow(train))
alfas_index<-mod_ksvm@alphaindex %>% unlist

alfa_zero[alfas_index]<-mod_ksvm@alpha %>% unlist
w_lol<-crossprod(alfa_zero,ifelse(train[,3]==1,1,-1)*as.matrix(train[,-3]))



library(e1071)
mod_e1071<-svm(classes~.,data=train,cost=1,kernel='linear',gamma = 1)

formula<-classes~.
kernel<-'linear'
scale=TRUE
threshold=0.001
gamma=1
C=1


svm<-function(  train,
                test,
                formula,
                kernel,
                scale=TRUE,
                threshold=0.001,
                gamma=1,
                C=1){
      
      #Defining the kernel functions
      polinomial_kernel<-function(X,Y,gamma,constant=0){
            poli_kern<- (tcrossprod(X,Y)+constant)^gamma
            return(poli_kern)
      }
      
      sigmoidal_kernel<-function(X,Y,gamma,constant=0){
            sig_kern<-tanh(gamma*tcrossprod(X,Y)+constant)
            return(sig_kern)
      }
      
      gaussian_kernel<-function(X,Y,gamma){
            dist_matrix<- as.matrix(dist(rbind(X,Y),method = 'euclidean'))[1:(nrow(X)),(nrow(train)+1):(nrow(X)+nrow(Y))]
            gaussian_kern<-exp(-gamma*dist_matrix^2)
            return(gaussian_kern)
      }
      
      
      exp_kernel<-function(X,Y,gamma){
            dist_matrix<- as.matrix(dist(rbind(X,Y),method = 'euclidean'))[1:(nrow(X)),(nrow(train)+1):(nrow(X)+nrow(Y))]
            exp_kern<-exp(-gamma*dist_matrix)
            return(exp_kern)
      }
      
      cauchy_kernel<-function(X,Y,gamma){
            dist_matrix<- as.matrix(dist(rbind(X,Y),method = 'euclidean'))[1:(nrow(X)),(nrow(train)+1):(nrow(X)+nrow(Y))]
            cauchy_kern<-1/(1+2*gamma*(dist_matrix^2))
            return(cauchy_kern)
      }
      
      laplacian_kernel<-function(X,Y,gamma){
            dist_matrix<- as.matrix(dist(rbind(X,Y),method = 'euclidean'))[1:(nrow(X)),(nrow(train)+1):(nrow(X)+nrow(Y))]
            laplacian_kern<-exp(-(dist_matrix)*sqrt(2*gamma))
            return(laplacian_kern)
      }
      
      
      noise<-1e-6
      var_names<-all.vars(formula)
      
      
      if(scale){
            #Scale Train
            mean_train<-apply(train[,colnames(train)!=var_names[1]],2,mean)
            sd<-apply(train[,colnames(train)!=var_names[1]],2,sd)
            train[,colnames(train)!=var_names[1]]<-t(apply(train[,colnames(train)!=var_names[1]],1,function(x){(x-mean_train)/sd}))
            
            #Scale test
            test[,colnames(test)!=var_names[1]]<-t(apply(test[,colnames(test)!=var_names[1]],1,function(x){(x-mean_train)/sd}))
            
            
      }
      #Setting the vars names
      
      X_train<-as.matrix(train[,colnames(train)!=var_names[1]])
      X_test<-as.matrix(test[,colnames(train)!=var_names[1]])
      
      
      n <- nrow(X_train)
      p <- ncol(X_train)+1
      n_test<-nrow(X_test)
      
      Y_train<-ifelse(train[,var_names[1]]==levels(train[,var_names[1]])[1],1,-1)
      Y_test<-test[,var_names[1]]
      
      
      #Here the vector b will be dimension 2n x 1, therefore, before that we will define an auxliar matrix, in way that the multiplication
      # of the vector b'A, gives the scalar on the left.
      A<-diag(Y_train,nrow=nrow(train),ncol=nrow(train)) # Auxiliar Diagonal Matrix Diagonal Matrix nxn of y
      
      

      if(kernel=="linear"){
            X_aux<-tcrossprod(X_train,X_train)
            X_aux_test<-tcrossprod(X_train,X_test)
            
      }else if(kernel=="polynomial"){
            X_aux<-polinomial_kernel(X=X_train,Y = X_train,gamma = gamma)
            X_aux_test<-polinomial_kernel(X=X_train,Y=X_test,gamma = gamma)
            
            
      }else if(kernel=="sigmoidal"){
            X_aux<-sigmoidal_kernel(X=X_train,Y = X_train,gamma = gamma)
            X_aux_test<-sigmoidal_kernel(X=X_train,Y=X_test,gamma = gamma)
            
            
      }else if(kernel=="gaussian"){
            X_aux<-gaussian_kernel(X=X_train,Y = X_train,gamma = gamma)
            X_aux_test<-gaussian_kernel(X=X_train,Y=X_test,gamma = gamma)
            
      }else if(kernel=="exponential"){
            X_aux<-exp_kernel(X=X_train,Y = X_train,gamma = gamma)
            X_aux_test<-exp_kernel(X=X_train,Y=X_test,gamma = gamma)
            
      }else if(kernel=="cauchy"){
            X_aux<-cauchy_kernel(X=X_train,Y = X_train,gamma = gamma)
            X_aux_test<-cauchy_kernel(X=X_train,Y=X_test,gamma = gamma)
            
      }else if(kernel=="laplacian"){
            X_aux<-laplacian_kernel(X=X_train,Y = X_train,gamma=gamma)
            X_aux_test<-laplacian_kernel(X=X_train,Y=X_test,gamma=gamma)
            
      }else{
            return("Invalid Kernel Function")
      }
      
      
      D<-A%*%tcrossprod(X_aux,A) #Dimension 2n x 2n, once that 2n variables will be optimized
      
      d<-matrix(rep(1,n),ncol=1)#Matriz column that need to multiply the vector d, has to have dimension 2n
      
      #Agora vamos escrever a matriz das
      B_1<-matrix(Y_train,nrow=1)#First restrictio of the sum
      
      
      B_4<-diag(1,n)#Restriction of alpha being greather than zero
      # B_4<-cbind(B_4,matrix(0,nrow=n,ncol=n))
      
      B_3<-diag(-1,n)#Restriction of alpha being less than C 
      # B_3<-cbind(B_3,matrix(0,nrow=n,ncol=n))
      
      
      B<-rbind(B_1,B_4,B_3)
      
      b_0<-c(0,rep(0,n),rep(-C,n)) # Right side of restrictions
      
      B<-t(B)
      
      library(quadprog)
      
      #Aqui ? o resultado final
      res<-solve.QP(D+diag(noise,nrow=nrow(D)),d,B,b_0,meq = 1)$solution #Here meq means that the first restriction is a equalitty.
      # res<-kernlab::ipop(c=c(Y_train-e,-Y_train-e),
      #                    H=D,A=c(rep(1,n),rep(-1,n)),
      #                    b=0,
      #                    l=rep(0,2*nrow(train)),u=rep(C,2*nrow(train)),r=0,sigf = 7)@primal
      
      alfa_total<-res
      
      #Probabaly notnecessary
      # alfa_up<-which(alfa_total>threshold)
      # alfa_down<-which(alfa_total< -threshold)
      
      alfa_total_aux<-which(abs(alfa_total)>=threshold)
      alfa_sv<-rep(0,length(alfa_total),nrow=nrow(train))
      alfa_sv[alfa_total_aux]<-alfa_total[alfa_total_aux]
      
      
      #A PORRA DE UM MENOS
      # alfa_sv<- -alfa_sv
      w<-crossprod(alfa_sv,Y_train*X_train)
      
      
      #Calculating b
      
      b_aux<-which(alfa_sv>0 & alfa_sv<C)
      arg<-which.max(alfa_sv)
      
      b<- (((1-(Y_train*tcrossprod(X_train,w)))/Y_train)) %>% median
      b
      

      
      #calculating Predicted Values for test and traininig
      y_pred_train<-alfa_sv%*%X_aux %>% as.vector
      y_pred_train<-sign(y_pred_train+b)
      
      
      y_pred_test<-alfa_sv%*%X_aux_test %>% as.vector
      y_pred_test<-sign(y_pred_test+b)
           #========================================================
      
      return(list(predict_train=y_pred_train,predict_test=y_pred_test))
}