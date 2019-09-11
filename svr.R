#Setting training and and test toy example
set.seed(42)
split_ratio<-0.7
data<-cars
train_index<-sample(1:nrow(data),round(split_ratio*nrow(data)))
train<-data[train_index,]
test<-data[-train_index,]



svr<-function(  train,
                test,
                formula,
                kernel,
                epsilon=0.1,
                scale=TRUE,
                threshold=0.001,
                gamma=1,
                C=1){

      noise<-1e-4
      
      if(scale){
            #Scale Train
            mean_train<-apply(train,2,mean)
            sd<-apply(train,2,sd)
            train<-t(apply(train,1,function(x){(x-mean_train)/sd}))

            #Scale test
            test<-t(apply(test,1,function(x){(x-mean_train)/sd}))
            
            
      }
      #Setting the vars names
      var_names<-all.vars(formula)
      
      X_train<-matrix(train[,colnames(train)!=var_names[1]])
      X_test<-matrix(test[,colnames(train)!=var_names[1]])
      
      
      n <- nrow(X_train)
      p <- ncol(X_train)
      n_test<-nrow(X_test)
      
      Y_train<-train[,var_names[1]]
      Y_test<-test[,var_names[1]]
      
      
      #Here the vector b will be dimension 2n x 1, therefore, before that we will define an auxliar matrix, in way that the multiplication
      # of the vector b'A, gives the scalar on the left.
      A_1<-diag(1,nrow=nrow(train),ncol=nrow(train)) # Diagonal Matrix nxn of 1
      A_2<-diag(-1,nrow=nrow(train))#Diagonal Matriz nxn of -1
      
      A<-rbind(A_1,A_2)#Auxiliar Matrix 2nxn, multipluing by W gives the first difference
      
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
      
      d<-matrix(c(Y_train-epsilon,-Y_train-epsilon),ncol=1)#Matriz column that need to multiply the vector d, has to have dimension 2n
      
      #Agora vamos escrever a matriz das
      B_1<-matrix(c(rep(1,n),rep(-1,n)),nrow=1)#First restrictio of the sum
      
      
      B_4<-diag(1,n)#Restriction of a+ being greather than zero
      B_4<-cbind(B_4,matrix(0,nrow=n,ncol=n))
      
      B_6<-diag(1,n)#Restriction of a- being greather than
      B_6<-cbind(matrix(0,nrow=n,ncol=n),B_6) 
      
      B_3<-diag(-1,n)#Restriction of a+ being less than C 
      B_3<-cbind(B_3,matrix(0,nrow=n,ncol=n))
      
      B_5<-diag(-1,n)#Restriction of a- being less than C 
      B_5<-cbind(matrix(0,nrow=n,ncol=n),B_5) 
      
      B<-rbind(B_1,B_4,B_6,B_3,B_5)
      
      b_0<-c(0,rep(0,n),rep(0,n),rep(-C,n),rep(-C,n)) #Lado direito das restri??es
      
      B<-t(B)
      
      library(quadprog)
      
      #Aqui ? o resultado final
      res<-solve.QP(D+diag(noise,nrow=nrow(D)),d,B,b_0,meq = 1)$solution #Here meq means that the first restriction is a equalitty.
      # res<-kernlab::ipop(c=c(Y_train-e,-Y_train-e),
      #                    H=D,A=c(rep(1,n),rep(-1,n)),
      #                    b=0,
      #                    l=rep(0,2*nrow(train)),u=rep(C,2*nrow(train)),r=0,sigf = 7)@primal
      
      alfa<-res[1:n]
      alfa_star<-res[(nrow(train)+1):(2*nrow(train))]
      
      
      alfa_total<-alfa_star-alfa
      
      #Probabaly notnecessary
      # alfa_up<-which(alfa_total>threshold)
      # alfa_down<-which(alfa_total< -threshold)
       
      alfa_total_aux<-which(abs(alfa_total)>=threshold)
      alfa_sv<-rep(0,length(alfa_total),nrow=nrow(train))
      alfa_sv[alfa_total_aux]<-alfa_total[alfa_total_aux]
      
      
      #A PORRA DE UM MENOS
      alfa_sv<- -alfa_sv
      w<-alfa_sv%*%X_train
      
      
      b1<-Y_train-alfa_sv%*%X_aux+epsilon
      b2<-Y_train-alfa_sv%*%X_aux-epsilon
      
      
      
      # b1<-Y_train[alfa_up]-w%*%X_train[alfa_up,]+e
      # b2<-Y_train[alfa_down]-w%*%X_train[alfa_down,]-e
      # 
      b<-median(c(b1,b2))
      b
      
      #calculating Predicted Values for test and traininig
      y_pred_train<-alfa_sv%*%X_aux %>% as.vector
      y_pred_train<-y_pred_train+b
      y_pred_train<-as.vector(y_pred_train)*sd[var_names[1]]+mean_train[var_names[1]]
      
      
      y_pred_test<-alfa_sv%*%X_aux_test %>% as.vector
      y_pred_test<-y_pred_test+b
      y_pred_test<-as.vector(y_pred_test)*sd[var_names[1]]+mean_train[var_names[1]]
      
      #========================================================
      
      return(list(predict_train=y_pred_train,predict_test=y_pred_test))
}


my_svr<-svr(train = train,test = test,formula = dist~.,kernel = 'linear',epsilon = 0.1,
            scale = TRUE,threshold = 0.001,gamma = 1,C = 1)


#Comparing with kernlab
library(kernlab)
mod_ksvm<-kernlab::ksvm(dist~.,data=train,C=1,epsilon=0.1,kernel='vanilladot')
pred_ksvm<-predict(mod_ksvm,newdata=train) %>% as.vector()

(pred_ksvm-my_svr$predict_train) %>% mean


#Comparing with e1071
library(e1071)
mod_e1071<-svm(dist~.,data=train,cost=1,epsilon=0.1,kernel='linear')
pred_e171<-predict(mod_e1071,newdata=test)

(pred_e171-my_svr$predict_test) %>% mean
